//file:noinspection GrMethodMayBeStatic
//file:noinspection GroovyUnusedAssignment
//file:noinspection GrDeprecatedAPIUsage

import groovy.transform.ImmutableOptions
import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.Geometry
import org.locationtech.jts.geom.GeometryFactory
import org.locationtech.jts.geom.LineString
import org.locationtech.jts.geom.Point
import org.locationtech.jts.linearref.LengthIndexedLine
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.interfaces.ROI
import qupath.lib.roi.GeometryTools

import static qupath.lib.scripting.QP.addObjects
import static qupath.lib.scripting.QP.getAnnotationObjects
import static qupath.lib.scripting.QP.getCurrentImageData
import static qupath.lib.scripting.QP.makeRGB
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
    def singleTissueSlide = tissuesWithLine.size() == 1
    for (TissueWithLines tissueAndLines : tissuesWithLine) {
        try {
            createGradientAnnotations(tissueAndLines, singleTissueSlide ? '' : "${coreIndex}_", liverBands, tumorBands)
        } catch (LocalRuntimeException e) {
            print("ERROR: coreIndex [${coreIndex}] for tissue [${tissueAndLines.tissue().getID()}] failed with error:\n${e.message}")
        }
        coreIndex++
    }

}

def createGradientAnnotations(TissueWithLines tissueAndLines, String annotationPrefix, List<Integer> liverBands, List<Integer> tumorBands) {
    List<PathObject> tissueAnnotations = []
    def (ROI liver, ROI tumor, ROI capsule) = getSeparatedTissueParts(tissueAndLines)
    tissueAnnotations << getAnnotation(tumor, "${annotationPrefix}tumor", makeRGB(150, 150, 0))
    tissueAnnotations << getAnnotation(liver, "${annotationPrefix}liver", makeRGB(150, 150, 0))

    def liverWC = unionROI(liver, capsule)
    def tumorWC = unionROI(tumor, capsule)
    if (capsule != null) {
        tissueAnnotations << getAnnotation(capsule, "${annotationPrefix}capsule_whole", makeRGB(0, 150, 0))
        def (midLine, liverCapsule, tumorCapsule) = splitCapsuleInHalves(capsule, tissueAndLines.lines(), tumor)
        tissueAnnotations << getAnnotation(midLine, "${annotationPrefix}capsule_midline", makeRGB(0, 150, 0))
        tissueAnnotations << getAnnotation(tumorCapsule, "${annotationPrefix}tumor_capsule", makeRGB(0, 150, 0))
        tissueAnnotations << getAnnotation(liverCapsule, "${annotationPrefix}liver_capsule", makeRGB(0, 150, 0))
    }
    def (biggestLiverExpansion, liverExpansionAnnotations) = annotateHalfWithExpansions("${annotationPrefix}tumor", liverWC, tumorWC, liverBands)
    tissueAnnotations.addAll(liverExpansionAnnotations)
    tissueAnnotations << getAnnotation(createCentralROI(tumor, biggestLiverExpansion), "${annotationPrefix}tumor_central", makeRGB(0, 150, 0))

    def (biggestTumorExpansion, tumorExpansionAnnotations) = annotateHalfWithExpansions("${annotationPrefix}liver", tumorWC, liverWC, tumorBands)
    tissueAnnotations.addAll(tumorExpansionAnnotations)
    tissueAnnotations << getAnnotation(createCentralROI(liver, biggestTumorExpansion), "${annotationPrefix}liver_central", makeRGB(0, 150, 0))
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

PathObject getAnnotation(ROI roi, String name = 'unnamed annotation', Integer color = makeRGB(255, 255, 0)) {
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
        annotateIgnoredGeometries(ignoredGeometries)

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
    annotateIgnoredGeometries(ignoredGeometries)
    return [liver, tumor, capsule].collect { toRoi(it) }
}

void annotateIgnoredGeometries(List<Geometry> ignoredGeometries) {
    def annotations = ignoredGeometries.collect { getAnnotation(toRoi(it), "00_ignored", makeRGB(200, 0, 100)) }
    addObjects(annotations.findAll { it.ROI.isLine() || it.ROI.getArea() > 0 })
}

def ALLOW_NO_CALIBRATION() { return true }

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
        midline = createMidlineStringV4(liverLine.ROI.geometry as LineString, tumorLine.ROI.geometry as LineString, capsule.geometry)
        if (midline == null) {
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

Geometry createMidlineStringV4(LineString line1, LineString line2, Geometry capsuleGeometry) {
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
        print('WARN: reversing line 2!')
        line2StartPoint = line2EndPoint
        line2 = line2.reverse()
    }

    def linesStartXDiff = line2StartPoint.getX() - line1StartPoint.getX()
    def linesStartYDiff = line2StartPoint.getY() - line1StartPoint.getY()
    def xCoefficient = linesStartYDiff / (Math.abs(linesStartYDiff) + Math.abs(linesStartXDiff))
    def yCoefficient = linesStartXDiff / (Math.abs(linesStartYDiff) + Math.abs(linesStartXDiff))
    def MAX_STEP_SIZE_POWER = 3
    def START_STEP_SIZE_POWER = -1
    for (int stepSizePower = START_STEP_SIZE_POWER; stepSizePower <= MAX_STEP_SIZE_POWER; stepSizePower++) {
        def midLinePoints = calcMidline(line1, line2, capsuleGeometry, stepSizePower)
        def midLineString = lineFromPoints(midLinePoints, geomFactory)
        if (!midLineString.intersects(line1) && !midLineString.intersects(line2)) {
            return midLineString
        }
    }
    return null
}

List<Point> calcMidline(LineString line1, LineString line2, Geometry capsule, int stepSizePower) {

    def liLine1 = new LengthIndexedLine(line1)
    def liLine2 = new LengthIndexedLine(line2)
    def geomFactory = capsule.factory
    def MAX_POINTS = line1.numPoints * (2**(stepSizePower))
    List<Point> points = []
    points.push(midPoint(line1.getStartPoint(), line2.getStartPoint(), line1.factory))
    def stepSize = line1.length / MAX_POINTS
    int i = 0
    def within = false
    def line1Offset = 0
    def line2Offset = 0
    while (true) {
        if (i > MAX_POINTS) {
            print('START_STEP_SIZE: ' + stepSize)
            print('MAX_START_POINTS: ' + MAX_POINTS)
            return points
        }
        i += 1
        def locationInLine1 = stepSize * i + line1Offset
        def p1 = liLine1.extractPoint(locationInLine1)
        def locationInLine2 = stepSize * i + line2Offset
        def p2 = liLine2.extractPoint(locationInLine2)
        def crossLine = lineFromCoords([p1, p2] as Coordinate[], geomFactory)
        def line1Intersections = crossLine.intersection(line1)
        if (line1Intersections.numPoints > 1) {
            // take point with closest distance to other side
            def closestP1 = line1Intersections.getCoordinates().sort { it.distance(p2) }.first()
            crossLine = lineFromCoords([closestP1, crossLine.getCoordinateN(1)] as Coordinate[], geomFactory)
        }
        def line2Intersections = crossLine.intersection(line2)

        if (line2Intersections.numPoints > 1) {
            // take point with closest distance to other side
            def closestP2 = line2Intersections.getCoordinates().sort { it.distance(p1) }.first()
            crossLine = lineFromCoords([crossLine.getCoordinateN(0), closestP2] as Coordinate[], geomFactory)
        }
//        addObject(getAnnotation(toRoi(crossLine), 'line ' + i))
        def newMidPoint = crossLine.centroid

        line1Offset += getOffsetToClosestNextPoint(liLine1, locationInLine1, newMidPoint)
        line2Offset += getOffsetToClosestNextPoint(liLine2, locationInLine2, newMidPoint)
        points.push(newMidPoint)
        def newWithin = newMidPoint.within(capsule)
        def breakingOut = within && !newWithin
        if (breakingOut) {
            break
        }
        within = newWithin
    }
    points.push(midPoint(line1.getEndPoint(), line2.getEndPoint(), line1.factory))
    return points
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

ROI toRoi(Geometry geometry, ImagePlane imagePlane = getDefaultImagePlane()) {
    return GeometryTools.geometryToROI(geometry, imagePlane)
}

ImagePlane getDefaultImagePlane() {
    return getAnnotationObjects().find { it.ROI != null }.ROI.imagePlane
}

Point midPoint(Point p1, Point p2, GeometryFactory geomFactory) {
    return midPoint(p1.coordinate, p2.coordinate, geomFactory)
}

Point midPoint(Coordinate p1, Coordinate p2, GeometryFactory geomFactory) {
    return geomFactory.createLineString([p1, p2] as Coordinate[]).centroid
}

LineString lineFromPoints(List<Point> points, GeometryFactory geomFactory) {
    return geomFactory.createLineString(points.collect { it.coordinate } as Coordinate[])
}

LineString lineFromCoords(Coordinate[] coords, GeometryFactory geomFactory) {
    return geomFactory.createLineString(coords)
}

int getOffsetToClosestNextPoint(LengthIndexedLine lengthIndexedLine, double lengthIndex, Point midPoint, double stepSize = 1d) {
    def bestDist = Double.POSITIVE_INFINITY
    def MAX_STEPS = (int) lengthIndexedLine.endIndex / stepSize
    for (int i = 0; i < MAX_STEPS; i++) {
        def offset = i * stepSize
        def newDist = lengthIndexedLine.extractPoint(lengthIndex + offset).distance(midPoint.coordinate)
        if (newDist >= bestDist) {
            return (i - 1) * stepSize
        }
        bestDist = newDist
    }
    return 0
}