//file:noinspection GrMethodMayBeStatic

import groovy.transform.ImmutableOptions
import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.Geometry
import org.locationtech.jts.geom.GeometryFactory
import org.locationtech.jts.geom.LineString
import qupath.lib.common.ColorTools
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.interfaces.ROI
import qupath.lib.roi.GeometryTools

import static qupath.lib.scripting.QP.addObject
import static qupath.lib.scripting.QP.addObjects
import static qupath.lib.scripting.QP.getAnnotationObjects
import static qupath.lib.scripting.QP.getCurrentImageData
import static qupath.lib.scripting.QP.removeObjects

print('START: main')
main()
print('END: main')

class LocalRuntimeException extends RuntimeException {
    LocalRuntimeException(String message) {
        super(message)
    }
}

def REJECT_DIRTY_ANNOTATIONS() { return false }

def main() {
    List<Integer> liverBands = [100] * 5 + [500] * 5
    List<Integer> tumorBands = [100] * 5 + [500] * 5
    cleanupAutoAnnotations()
    List<TissueWithLines> tissuesWithLine = findTissueWithTumorLines()
    Integer coreIndex = 0
    for (TissueWithLines tissueAndLines : tissuesWithLine) {
        try {
            createGradientAnnotations(tissueAndLines, coreIndex, liverBands, tumorBands)
        } catch (LocalRuntimeException e) {
            print("ERROR: coreIndex [${coreIndex}] for tissue [${tissueAndLines.tissue().getID()}] failed with error:\n${e.message}")
        }
        coreIndex++
    }

}

def DEBUG_MODE() { return false }

def DEBUG_MODE_CAPSULE_DIRECTIONS() { return false }

def VISUALIZE_PATH_FINDING() { return true }

def createGradientAnnotations(TissueWithLines tissueAndLines, int coreIndex, List<Integer> liverBands, List<Integer> tumorBands) {
    List<PathObject> tissueAnnotations = []
    def (liver, tumor, capsule) = getSeparatedTissueParts(tissueAndLines)
    tissueAnnotations << getAnnotation(tumor, "${coreIndex}_tumor", makeRGB(150, 150, 0))
    tissueAnnotations << getAnnotation(liver, "${coreIndex}_liver", makeRGB(150, 150, 0))

    def liverWC = unionROI(liver, capsule)
    def tumorWC = unionROI(tumor, capsule)
    if (capsule != null) {
        tissueAnnotations << getAnnotation(capsule, "${coreIndex}_capsule_whole", makeRGB(0, 150, 0))
        def (midLine, liverCapsule, tumorCapsule) = splitCapsuleInHalves(capsule, tissueAndLines.lines(), tumor)
        tissueAnnotations << getAnnotation(midLine, "${coreIndex}_capsule_midline", makeRGB(0, 150, 0))
        tissueAnnotations << getAnnotation(tumorCapsule, "${coreIndex}_capsule_tumor", makeRGB(0, 150, 0))
        tissueAnnotations << getAnnotation(liverCapsule, "${coreIndex}_capsule_liver", makeRGB(0, 150, 0))
    }
    def (biggestLiverExpansion, liverExpansionAnnotations) = annotateHalfWithExpansions("${coreIndex}_tumor", liverWC, tumorWC, liverBands)
    tissueAnnotations.addAll(liverExpansionAnnotations)
    tissueAnnotations << getAnnotation(createCentralROI(tumor, biggestLiverExpansion), "${coreIndex}_tumor_central", makeRGB(0, 150, 0))

    def (biggestTumorExpansion, tumorExpansionAnnotations) = annotateHalfWithExpansions("${coreIndex}_liver", tumorWC, liverWC, tumorBands)
    tissueAnnotations.addAll(tumorExpansionAnnotations)
    tissueAnnotations << getAnnotation(createCentralROI(liver, biggestTumorExpansion), "${coreIndex}_liver_central", makeRGB(0, 150, 0))
    addObjects(tissueAnnotations.findAll { it.ROI.isLine() || it.ROI.getArea() > 0 })
}


def annotateHalfWithExpansions(String name, ROI startPolygon, ROI intersectingPolygon, List<Integer> bandsInMicrons) {
    def prevExpansion = startPolygon.geometry
    List<PathObject> newAnnotations = []
    def totalMicrons = 0
    for (int i = 1; i <= bandsInMicrons.size(); i++) {
        totalMicrons += bandsInMicrons[i - 1]
        def totalDistance = getDistance(totalMicrons)
        def expansion = startPolygon.geometry.buffer(totalDistance)
        def expansionBand = expansion.difference(prevExpansion)
        prevExpansion = expansion
        def intersection = expansionBand.intersection(intersectingPolygon.geometry)

        def bandColor = makeRGB(20 * i, 40 * i, 200 - 30 * i)
        newAnnotations << getAnnotation(toRoi(intersection), "${name}_${String.format('%04d', totalMicrons)}µm", bandColor)
    }
    return [prevExpansion, newAnnotations]
}

List<TissueWithLines> findTissueWithTumorLines() {
    Collection<PathObject> annotations = getAnnotationObjects()
    // Find required annotations
    return annotations.collect { getTissueWithLine(it) }.findAll { it !== null }
}

@ImmutableOptions(knownImmutableClasses = [PathObject])
record TissueWithLines(PathObject tissue, Collection<PathObject> lines) {}

TissueWithLines getTissueWithLine(PathObject candidateTissue) {
    assert candidateTissue != null
    if (!candidateTissue.ROI.isArea() || candidateTissue.classifications.contains('Auto')) {
        return null
    }

    print "Found tissue annotation: ${candidateTissue?.name}"

    Collection<PathObject> annotations = getAnnotationObjects()
    Collection<PathObject> tissueLines = findLinesForTissue(annotations, candidateTissue)
    if (!(tissueLines?.size() > 0)) {
        return null
    }
    return new TissueWithLines(candidateTissue, tissueLines)
}

Collection<PathObject> findLinesForTissue(Collection<PathObject> annotations, PathObject tissue) {
    Collection<PathObject> lines = annotations.findAll { it.ROI.isLine() && it.ROI.geometry.intersects(tissue.ROI.geometry) }
    if (lines.size() == 0) {
        return null
    } else if (lines.size() == 1) {
        print "Found tumor line for tissue: [${tissue?.getID()}]"
    } else if (lines.size() == 2) {
        print "Found tumor and liver line for tissue: [${tissue?.getID()}]"
    } else {
        throw new LocalRuntimeException("More than 2 line annotations found for tissue [${tissue?.getID()}]")
    }
    return lines
}

PathObject getAnnotation(ROI roi, String name, Integer color = makeRGB(255, 255, 0)) {
    def newAnnotation = PathObjects.createAnnotationObject(roi)
    newAnnotation.setClassifications(['Auto'])
    newAnnotation.setColor(color)
    newAnnotation.setName(name)
    Collection<PathObject> annotations = getAnnotationObjects()
    def existing = annotations.find { it.name == name }
    if (existing != null) {
        print "Removing existing annotation [$name]"
        removeObject(existing, false)
    }
    print "Successfully created ROI [$name] annotation"
    return newAnnotation
}

Tuple<ROI> getSeparatedTissueParts(TissueWithLines tissueAndLines) {
    def (tissue, tissueSplitLines) = [tissueAndLines.tissue(), tissueAndLines.lines()]
    def halvesGeometries = GeometryTools.splitGeometryByLineStrings(tissue.ROI.geometry, tissueSplitLines.collect { it.ROI.geometry })
    def tumors = halvesGeometries.findAll { isTumor(it) }.sort { it.area }

    if (tumors.size() === 0) {
        addObjects(halvesGeometries.collect { getAnnotation(toRoi(it), "00_debug_tissue_without_tumor_halves", makeRGB(255, 50, 50)) })
        throw new LocalRuntimeException("No `Tumor` annotation found within the bounds of tissue [${tissue?.getID()}], please add an annotation with `Tumor` classification manually within one the tissue somewhere.")
    }
    def tumor = tumors.last
    if (tissueSplitLines.size() == 1) {
        if (halvesGeometries.size() !== 2) {
            if (halvesGeometries.size() < 2 || REJECT_DIRTY_ANNOTATIONS()) {
                addObjects(halvesGeometries.collect { getAnnotation(toRoi(it), "00_debug", makeRGB(255, 50, 50)) })
                throw new LocalRuntimeException("Expected 2 halves, but got ${halvesGeometries.size()}, please check the 00_debug annotation, probably the smallest one is the issue")
            } else {
                print("WARNING: Expected 2 halves, but got ${halvesGeometries.size()}, this could give weird results, please check the annotations, setting REJECT_DIRTY_ANNOTATIONS to true can help with that.")
            }
        }
        def liver = halvesGeometries.findAll() { !isTumor(it) }.sort { it.area }.last

        def ignoredGeometries = halvesGeometries.findAll { ![liver, tumor].contains(it) }
        annotateIgnoredGeometries(ignoredGeometries, tissue)

        return [liver, tumor].collect { toRoi(it) }
    }

    // We have more than 1 line, we need to find the capsule
    if (halvesGeometries.size() !== 3) {
        if (halvesGeometries.size() < 3 || REJECT_DIRTY_ANNOTATIONS()) {
            addObjects(halvesGeometries.collect { getAnnotation(toRoi(it), "00_debug", makeRGB(255, 50, 50)) })
            throw new LocalRuntimeException("Expected 3 halves, but got ${halvesGeometries.size()}")
        } else {
            print("WARNING: Expected 3 halves, but got ${halvesGeometries.size()}, this could give weird results, please check the annotations, setting REJECT_DIRTY_ANNOTATIONS to true can help with that.")
        }
    }
    def capsuleGeometries = halvesGeometries.findAll { it.touches(tumor) }
    if (capsuleGeometries.size() < 1) {
        addObjects([tissue] + tissueSplitLines.collect { getAnnotation(it.ROI, "00_debug", makeRGB(255, 50, 50)) })
        throw new LocalRuntimeException("Expected 1 or more geometry that touches the tumor geometry, but got ${capsuleGeometries.size()}")
    }
    def capsule = mergeGeometries(capsuleGeometries)

    def liver = halvesGeometries.findAll { !isTumor(it) && !capsuleGeometries.contains(it) }.sort { it.area }.last
    def ignoredGeometries = halvesGeometries.findAll { !([liver, tumor] + capsuleGeometries).contains(it) }
    annotateIgnoredGeometries(ignoredGeometries, tissue)
    return [liver, tumor, capsule].collect { toRoi(it) }
}

void annotateIgnoredGeometries(List<Geometry> ignoredGeometries, tissue) {
    def annotations = ignoredGeometries.collect { getAnnotation(toRoi(it), "00_ignored", makeRGB(200, 0, 100)) }
    addObjects(annotations.findAll { it.ROI.isLine() || it.ROI.getArea() > 0 })
}

def ALLOW_NO_CALIBRATION() { return false }

double getDistance(double microns) {
    def imageData = getCurrentImageData()
    def server = imageData.getServer()
// We need the pixel size
    def cal = server.getPixelCalibration()
    if (!cal.hasPixelSizeMicrons()) {

        if (cal.getPixelWidthUnit() == 'px' && ALLOW_NO_CALIBRATION()) {
            print 'Warning!!! Going through with pixels instead of microns for debugging purposes'
            return microns / cal.getAveragedPixelSize()
        } else {
            print 'We need the pixel size information here!'
            print 'cal.getPixelWidth() ' + cal.getPixelWidth()
            print 'cal.getPixelWidthMicrons() ' + cal.getPixelWidthMicrons()
            print 'cal.getPixelWidthUnit(): ' + cal.getPixelWidthUnit()
            print 'cal.getAveragedPixelSizeMicrons(): ' + cal.getAveragedPixelSizeMicrons()
        }
        throw new LocalRuntimeException("We need the pixel size information here!")
    }
    return microns / cal.getAveragedPixelSizeMicrons()
}

boolean isTumor(Geometry geom) {
    Collection<PathObject> annotations = getAnnotationObjects()
    return annotations.find { it.classifications.contains('Tumor') && geom.intersects(it.ROI.geometry) } != null
}

ROI createCentralROI(ROI startRoi, Geometry biggestExpansion) {
    def centralGeometry = startRoi.geometry.difference(biggestExpansion)
    return toRoi(centralGeometry)
}

void cleanupAutoAnnotations() {
    Collection<PathObject> annotations = getAnnotationObjects()
    removeObjects(annotations.findAll { it.classifications.contains('Auto') }, false)
}


ROI unionROI(ROI roi1, ROI roi2) {
    if (roi2 == null) {
        return roi1
    }
    if (roi1 == null) {
        return roi2
    }
    return toRoi(roi1.geometry.union(roi2.geometry))
}

Tuple<ROI> splitCapsuleInHalves(ROI capsule, Collection<PathObject> lines, ROI tumor) {
    def capsuleGeometry = capsule.geometry
    def tumorLine = lines[0].ROI.geometry.touches(tumor.geometry) ? lines[0] : lines[1]
    def liverLine = lines[0] == tumorLine ? lines[1] : lines[0]
    Geometry midline = null
    try {
        if (tumorLine.ROI.geometry.intersects(liverLine.ROI.geometry)) {
            throw new LocalRuntimeException('Liver line is intersecting tumor line')
        }
        midline = createMidlineStringV4(liverLine.ROI.geometry, tumorLine.ROI.geometry, capsule.geometry)
        if (midline == null || midline.intersects(tumorLine.ROI.geometry) || midline.intersects(liverLine.ROI.geometry)) {
            throw new LocalRuntimeException('Could not find a non-intersecting midline for the capsule')
        }
        if (!midline.intersects(capsuleGeometry)) {
            throw new LocalRuntimeException('Expected (midline.intersects(capsuleGeometry))')
        }
        def capsuleParts = GeometryTools.splitGeometryByLineStrings(capsuleGeometry, [midline])
        if (capsuleParts.size() < 2) {
            throw new LocalRuntimeException('Expected (capsuleParts.size() >= 2)')
        }
        def tumorCapsule = mergeGeometries(capsuleParts.findAll { it.touches(tumor.geometry) || it.intersects(tumor.geometry) })
        if (tumorCapsule == null) {
            addObjects(capsuleParts.collect { getAnnotation(toRoi(it), "00_debug_capsule_part", makeRGB(255, 50, 50)) })
            throw new LocalRuntimeException('Expected (tumorCapsule != null)')
        }
        def liverCapsule = mergeGeometries(capsuleParts.findAll { !it.touches(tumor.geometry) })
        if (liverCapsule == null) {
            throw new LocalRuntimeException('Expected (liverCapsule != null)')
        }
        return [midline, liverCapsule, tumorCapsule].collect { toRoi(it) }
    } catch (e) {
        addDebugAnnotations(tumorLine, liverLine, midline, capsule)
        throw e
    }
}

def addDebugAnnotations(PathObject tumorLine, PathObject liverLine, Geometry midline, ROI capsule) {
    List<PathObject> debugAnnotations = []
    debugAnnotations << getAnnotation(tumorLine.ROI, "00_debug_tumorLine", makeRGB(255, 50, 50))
    debugAnnotations << getAnnotation(liverLine.ROI, "00_debug_liverLine", makeRGB(255, 50, 50))
    if (midline !== null) {
        debugAnnotations << getAnnotation(toRoi(midline), "00_debug_midline", makeRGB(255, 50, 50))
    }
    debugAnnotations << getAnnotation(capsule, "00_debug_capsuleGeometry", makeRGB(255, 50, 50))
    addObjects(debugAnnotations)
}

Geometry createMidlineStringV4(Geometry line1, Geometry line2, Geometry capsuleGeometry) {
    if (!(line1 instanceof LineString)) {
        throw new IllegalArgumentException("line1 must be a LineString")
    }
    if (!(line2 instanceof LineString)) {
        throw new IllegalArgumentException("line2 must be a LineString")
    }

    def geomFactory = line1.getFactory()
    // First, we create straight line between the start and end of a line that is between the midpoints between the start and end of the two lines

    def line2StartPoint = line2.getStartPoint()
    def line2EndPoint = line2.getEndPoint()

    def line1StartPoint = line1.getStartPoint()
    def line1XCoefficient = line1.getPointN(1).getX() - line1.getPointN(0).getX()
    if (line2StartPoint.distance(line1StartPoint) > line2EndPoint.distance(line1StartPoint)) {
        line2StartPoint = line2EndPoint
    }
    def linesStartXDiff = line2StartPoint.getX() - line1StartPoint.getX()
    def linesStartYDiff = line2StartPoint.getY() - line1StartPoint.getY()
    def xCoefficient = linesStartYDiff / (Math.abs(linesStartYDiff) + Math.abs(linesStartXDiff))
    def yCoefficient = linesStartXDiff / (Math.abs(linesStartYDiff) + Math.abs(linesStartXDiff))
    if (Math.abs(xCoefficient - line1XCoefficient) > Math.abs(xCoefficient + line1XCoefficient)) {
        // Select the direction sense that matches up most with the fist segment of line1
        xCoefficient = -xCoefficient
        yCoefficient = -yCoefficient
    }
    // Then we create a LengthIndexedLine from the reference line
    def midLineStart = new Coordinate((line1StartPoint.x + line2StartPoint.x) / 2.0, (line1StartPoint.y + line2StartPoint.y) / 2.0)
    def MAX_POWER = 6
    Geometry midLine = null
    def annotations = []
    def refLinePoints = Math.max(line1.numPoints, line2.numPoints)
    List<PathObject> segmentAnnotations = []
    Collection<PathObject> renderedSegments = []
    for (int midLineResolutionPower = 0; midLineResolutionPower <= MAX_POWER; midLineResolutionPower++) {
        def sampleCount = refLinePoints * (2**(midLineResolutionPower - 3))
        def partSize = 2 * line1.getLength() / sampleCount
        print('Trying to create a capsule midline with resolution power ' + midLineResolutionPower + ', sampleCount ' + sampleCount)
        def midPoints = [midLineStart]
        xCoefficient = linesStartYDiff / (Math.abs(linesStartYDiff) + Math.abs(linesStartXDiff))
        yCoefficient = linesStartXDiff / (Math.abs(linesStartYDiff) + Math.abs(linesStartXDiff))
        // We loop over the sampleCount,
        for (int i = 0; i <= sampleCount; i++) {
            if ((i % 10) == 0) Thread.sleep(0) // Enable killing the process for if it takes too long

            def prev = midPoints.last
            def newPoint = new Coordinate(prev.getX() + xCoefficient * partSize, prev.getY() + yCoefficient * partSize)

            if (DEBUG_MODE()) annotations << getAnnotation(toRoi(geomFactory.createPoint(newPoint)), "00_debug_point_start_" + i, ColorTools.makeRGBA(20, 20, 20, 100))
            if (capsuleGeometry.contains(geomFactory.createPoint(prev))) {
                def newMidPoint = findNewMidPoint(prev, newPoint, line1, line2, geomFactory, annotations, i)

                def newPointIsOutside = !capsuleGeometry.contains(geomFactory.createPoint(newMidPoint))
                def insideToOutsideCapsule = capsuleGeometry.contains(geomFactory.createPoint(prev)) && newPointIsOutside
                if (VISUALIZE_PATH_FINDING()) segmentAnnotations << getAnnotation(toRoi(geomFactory.createLineString([prev, newMidPoint] as Coordinate[])),
                        "00_debug_seg_" + i, ColorTools.makeRGBA(i % 2 == 0 ? 255 : 0, i % 2 == 0 ? 255 : 0, i % 2 == 0 ? 255 : 0, 255)
                )


                if (newMidPoint !== null) {
//                print("Coefficients updated from [${xCoefficient}] [${yCoefficient}]")
                    double newXDiff = newMidPoint.getX() - prev.getX()
                    double newYDiff = newMidPoint.getY() - prev.getY()
                    xCoefficient = newXDiff / (Math.abs(newYDiff) + Math.abs(newXDiff))
                    yCoefficient = newYDiff / (Math.abs(newYDiff) + Math.abs(newXDiff))
//                print("Coefficients updated to [${xCoefficient}] [${yCoefficient}]")
                    midPoints << newMidPoint
                } else {
                    midPoints << newPoint
                }
                if (insideToOutsideCapsule) {
                    break
                }
            } else {
                midPoints << newPoint
            }
            if (i % 10 == 0) {
                if (VISUALIZE_PATH_FINDING()) {
                    addObjects(segmentAnnotations)
                    renderedSegments.addAll(segmentAnnotations)
                    segmentAnnotations.clear()
                }
                if (DEBUG_MODE()) {
                    addObjects(annotations)
                    annotations.clear()
                }
            }
        }
        if (renderedSegments.size() > 0) {
            removeObjects(renderedSegments, false)
        }
//        midPoints << refLineEnd
        if (midPoints.size() >= 2) {
            midLine = geomFactory.createLineString(midPoints as Coordinate[])
            if (!midLine.intersects(line1) && !midLine.intersects(line2)) {
                def capsuleCrossings = midLine.intersection(capsuleGeometry.boundary)
                if (capsuleCrossings.numPoints == 2) {
                    return midLine
                } else {
                    print("WARN midline finished but did not cross capsule in 2 points, but instead ${capsuleCrossings.numPoints}")
                }
            }
        }
        // If the midline intersects the lines, we need to try again with a higher resolution
    }
    return midLine
}

double calculateAngleDegrees(double opposite, double adjacent) {
    double radians = Math.atan(opposite / adjacent)
    return Math.toDegrees(radians)
}

Coordinate findNewMidPoint(Coordinate prev, Coordinate newPoint, LineString line1, LineString line2, GeometryFactory geomFactory, ArrayList<PathObject> annotations, int i) {
    def crossLineLength = line1.length
    double newXDiff = newPoint.getX() - prev.getX()
    double newYDiff = newPoint.getY() - prev.getY()
    if (newXDiff.isInfinite() || newYDiff.isInfinite() || newXDiff.isNaN() || newYDiff.isNaN()) {
        print('newPoint.getY()' + newPoint.getY() + ', i' + i)
        print('prev.getY()' + prev.getY() + ', i' + i)
        print('new newYDiff:' + newYDiff)
        throw new RuntimeException('Bad distances')
    }
    try {
        def prevAngle = (int) calculateAngleDegrees(-newYDiff, newXDiff)
        def betterPoint = findBestNewPoint(newPoint, crossLineLength, prevAngle, geomFactory, line1, line2, i)
//        annotations << getAnnotation(toRoi(geomFactory.createPoint(p1)), "00_debug_p1_" + i, makeRGB(255, 50, 50))
//        annotations << getAnnotation(toRoi(geomFactory.createPoint(p2)), "00_debug_p2_" + i, makeRGB(255, 50, 50))
        return betterPoint
    } catch (LocalRuntimeException e) {
        print('ERROR: error while trying to find new midpoint')
        print('ERROR: ' + e.message)
        if (DEBUG_MODE()) {
            throw e
        }
        print('WARN: ignoring error above')
    }
}

Coordinate findBestNewPoint(Coordinate newPoint, double orthLength, int prevAngle, GeometryFactory geomFactory, LineString line1, LineString line2, int i) {
    def minCrossLength = orthLength
    Tuple<Coordinate> bestPoints = [null, null]
    def ANGLE_STEP_SIZE = 20
    def MAX_ANGLE_DIFF = 100
    List<PathObject> debugAnnotations = []
    for (int dir = -1; dir <= 1; dir += 2) {
        for (int angleDiff = 0; angleDiff <= MAX_ANGLE_DIFF; angleDiff += ANGLE_STEP_SIZE) {
            angle = prevAngle + angleDiff * dir
            def xDir = Math.sin(angle * 2f * Math.PI / 360f)
            def yDir = Math.cos(angle * 2f * Math.PI / 360f)
            def firstOrthStart = new Coordinate(newPoint.x - xDir * orthLength, newPoint.y - yDir * orthLength)
            def firstOrthEnd = new Coordinate(newPoint.x + xDir * orthLength, newPoint.y + yDir * orthLength)

            LineString orthogonalLine = geomFactory.createLineString([firstOrthStart, firstOrthEnd] as Coordinate[])
            if (DEBUG_MODE_CAPSULE_DIRECTIONS()) debugAnnotations << getAnnotation(toRoi(orthogonalLine), "00_debug_direction_i${i}_a$angle", ColorTools.makeRGBA(20, 20, 20, 75))

            def p1 = selectClosestPoint(line1.intersection(orthogonalLine), newPoint)
            def p2 = selectClosestPoint(line2.intersection(orthogonalLine), newPoint)
            if (p1 != null && p2 != null) {
                def p1Dist = Math.abs(p1.distance(newPoint))
                def p2Dist = Math.abs(p2.distance(newPoint))
                def crossLength = p1Dist + p2Dist
                if (crossLength < minCrossLength) {
                    minCrossLength = crossLength
                    bestPoints = [p1, p2]
                }
            }
        }
    }
    if (DEBUG_MODE()) addObjects(debugAnnotations)
    def (p1, p2) = bestPoints
    if (p1 == null || p2 == null) {
        // This is the case of orthogonal line on the edges that is not intersecting one of the 2 lines, this is normal
        print("SKIPPING point $i, no intersection found p1 [$p1] p2 [$p2]")
        return null
    } else {
        LineString crossLine = geomFactory.createLineString([p1, p2] as Coordinate[])
//            annotations << getAnnotation(toRoi(crossLine), "00_debug_cross_" + i, ColorTools.makeRGBA(20, 20, 20, 75))
        if (DEBUG_MODE()) {
            addObject(getAnnotation(toRoi(crossLine), "00_debug_cross_" + i, ColorTools.makeRGBA(20, 60, 20, 75)))
        }   // Calculate the midpoint
        def newMidPoint = crossLine.centroid
//        annotations << getAnnotation(toRoi(newMidPoint), "00_debug_newMid" + i, makeRGB(255, 50, 50))
        return newMidPoint.coordinate
    }
}


Geometry mergeGeometries(List<Geometry> geometries) {
    if (geometries.size() == 0) {
        return null
    }
    def merged = geometries[0]
    for (int i = 1; i < geometries.size(); i++) {
        merged = merged.union(geometries[i])
    }
    return merged
}

Coordinate selectClosestPoint(Geometry geometry, Coordinate other) {
    def closestDist = Double.POSITIVE_INFINITY
    Coordinate closest = null
    if (geometry == null) {
        return null
    }

    def coordinates = geometry.getCoordinates()
    for (int i = 0; i < coordinates.length; i++) {
//        print('loop i ' + i)
        def option = coordinates[i]
        def newDist = Math.abs(option.distance(other))
//        print('loop newDist ' + newDist)
        if (closestDist > newDist) {
//            print('loop closestDist ' + newDist)
            closestDist = newDist
            closest = option
        }
    }
    return closest
}

def toRoi(Geometry geometry, imagePlane = getDefaultImagePlane()) {
    return GeometryTools.geometryToROI(geometry, imagePlane)
}

ImagePlane getDefaultImagePlane() {
    return getAnnotationObjects().find { it.ROI != null }.ROI.imagePlane
}

