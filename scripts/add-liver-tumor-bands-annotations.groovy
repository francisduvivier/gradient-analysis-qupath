//file:noinspection GrMethodMayBeStatic

import groovy.transform.ImmutableOptions
import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.Geometry
import org.locationtech.jts.geom.LineString
import org.locationtech.jts.linearref.LengthIndexedLine
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.roi.interfaces.ROI
import qupath.lib.roi.GeometryTools

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

def createGradientAnnotations(TissueWithLines tissueAndLines, int coreIndex, List<Integer> liverBands, List<Integer> tumorBands) {
    def (liver, tumor, capsule) = getSeparatedTissueParts(tissueAndLines)
    def liverWC = unionROI(liver, capsule)
    def tumorWC = unionROI(tumor, capsule)
    List<PathObject> tissueAnnotations = []
    if (capsule != null) {
        tissueAnnotations << getAnnotation(capsule, "${coreIndex}_capsule_whole", makeRGB(0, 150, 0))
        def (midLine, liverCapsule, tumorCapsule) = splitCapsuleInHalves(capsule, tissueAndLines.lines(), tumor)
        tissueAnnotations << getAnnotation(midLine, "${coreIndex}_capsule_midline", makeRGB(0, 150, 0))
        tissueAnnotations << getAnnotation(tumorCapsule, "${coreIndex}_capsule_tumor", makeRGB(0, 150, 0))
        tissueAnnotations << getAnnotation(liverCapsule, "${coreIndex}_capsule_liver", makeRGB(0, 150, 0))
    }
    tissueAnnotations << getAnnotation(liver, "${coreIndex}_liver", makeRGB(150, 150, 0))
    def (biggestLiverExpansion, liverExpansionAnnotations) = annotateHalfWithExpansions("${coreIndex}_tumor", liverWC, tumorWC, liverBands)
    tissueAnnotations.addAll(liverExpansionAnnotations)
    tissueAnnotations << getAnnotation(createCentralROI(tumor, biggestLiverExpansion), "${coreIndex}_tumor_central", makeRGB(0, 150, 0))

    tissueAnnotations << getAnnotation(tumor, "${coreIndex}_tumor", makeRGB(150, 150, 0))
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
        newAnnotations << getAnnotation(GeometryTools.geometryToROI(intersection, startPolygon.imagePlane), "${name}_${String.format('%04d', totalMicrons)}µm", bandColor)
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
    print "Successfully added ROI [$name] annotation"
    return newAnnotation
}

Tuple<ROI> getSeparatedTissueParts(TissueWithLines tissueAndLines) {
    def (tissue, tumorLines) = [tissueAndLines.tissue(), tissueAndLines.lines()]
    def halvesGeometries = GeometryTools.splitGeometryByLineStrings(tissue.ROI.geometry, tumorLines.collect { it.ROI.geometry })
    def tumors = halvesGeometries.findAll { isTumor(it) }.sort { it.area }

    if (tumors.size() === 0) {
        addObjects(halvesGeometries.collect { getAnnotation(GeometryTools.geometryToROI(it, tissue.ROI.imagePlane), "00_debug_tissue_without_tumor_halves", makeRGB(255, 50, 50)) })
        throw new LocalRuntimeException("No `Tumor` annotation found within the bounds of tissue [${tissue?.getID()}], please add an annotation with `Tumor` classification manually within one the tissue somewhere.")
    }
    def tumor = tumors.last
    if (tumorLines.size() == 1) {
        if (halvesGeometries.size() !== 2) {
            if (halvesGeometries.size() < 2 || REJECT_DIRTY_ANNOTATIONS()) {
                addObjects(halvesGeometries.collect { getAnnotation(GeometryTools.geometryToROI(it, tissue.ROI.imagePlane), "00_debug", makeRGB(255, 50, 50)) })
                throw new LocalRuntimeException("Expected 2 halves, but got ${halvesGeometries.size()}, please check the 00_debug annotation, probably the smallest one is the issue")
            } else {
                print("WARNING: Expected 2 halves, but got ${halvesGeometries.size()}, this could give weird results, please check the annotations, setting REJECT_DIRTY_ANNOTATIONS to true can help with that.")
            }
        }
        def liver = halvesGeometries.findAll() { !isTumor(it) }.sort { it.area }.last

        def ignoredGeometries = halvesGeometries.findAll { ![liver, tumor].contains(it) }
        annotateIgnoredGeometries(ignoredGeometries, tissue)

        return [liver, tumor].collect { GeometryTools.geometryToROI(it, tissue.ROI.imagePlane) }
    }

    // We have more than 1 line, we need to find the capsule
    if (halvesGeometries.size() !== 3) {
        if (halvesGeometries.size() < 3 || REJECT_DIRTY_ANNOTATIONS()) {
            addObjects(halvesGeometries.collect { getAnnotation(GeometryTools.geometryToROI(it, tissue.ROI.imagePlane), "00_debug", makeRGB(255, 50, 50)) })
            throw new LocalRuntimeException("Expected 3 halves, but got ${halvesGeometries.size()}")
        } else {
            print("WARNING: Expected 3 halves, but got ${halvesGeometries.size()}, this could give weird results, please check the annotations, setting REJECT_DIRTY_ANNOTATIONS to true can help with that.")
        }
    }
    def capsuleGeometries = halvesGeometries.findAll { it.touches(tumor) }
    if (capsuleGeometries.size() < 1) {
        addObjects([tissue] + tumorLines.collect { getAnnotation(it.ROI, "00_debug", makeRGB(255, 50, 50)) })
        throw new LocalRuntimeException("Expected 1 or more geometry that touches the tumor geometry, but got ${capsuleGeometries.size()}")
    }
    def capsule = mergeGeometries(capsuleGeometries)

    def liver = halvesGeometries.findAll { !isTumor(it) && !capsuleGeometries.contains(it) }.sort { it.area }.last
    def ignoredGeometries = halvesGeometries.findAll { !([liver, tumor] + capsuleGeometries).contains(it) }
    annotateIgnoredGeometries(ignoredGeometries, tissue)
    return [liver, tumor, capsule].collect { GeometryTools.geometryToROI(it, tissue.ROI.imagePlane) }
}

void annotateIgnoredGeometries(List<Geometry> ignoredGeometries, tissue) {
    def annotations = ignoredGeometries.collect { getAnnotation(GeometryTools.geometryToROI(it, tissue.ROI.imagePlane), "00_ignored", makeRGB(200, 0, 100)) }
    addObjects(annotations.findAll { it.ROI.isLine() || it.ROI.getArea() > 0 })
}

def ALLOW_NO_CALIBRATION() { return false }

double getDistance(double microns) {
    def imageData = getCurrentImageData()
    def server = imageData.getServer()
// We need the pixel size
    def cal = server.getPixelCalibration()
    if (!cal.hasPixelSizeMicrons()) {
        print 'We need the pixel size information here!'
        print 'cal.getPixelWidth() ' + cal.getPixelWidth()
        print 'cal.getPixelWidthMicrons() ' + cal.getPixelWidthMicrons()
        print 'cal.getPixelWidthUnit(): ' + cal.getPixelWidthUnit()
        print 'cal.getAveragedPixelSizeMicrons(): ' + cal.getAveragedPixelSizeMicrons()
        if (cal.getPixelWidthUnit() == 'px' && ALLOW_NO_CALIBRATION()) {
            print 'Warning!!! Going through with pixels instead of microns for debugging purposes'
            return microns / cal.getAveragedPixelSize()
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
    return GeometryTools.geometryToROI(centralGeometry, startRoi.imagePlane)
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
    return GeometryTools.geometryToROI(roi1.geometry.union(roi2.geometry), roi1.imagePlane)
}

Tuple<ROI> splitCapsuleInHalves(ROI capsule, Collection<PathObject> lines, ROI tumor) {
    def capsuleGeometry = capsule.geometry
    def tumorLine = lines[0].ROI.geometry.touches(tumor.geometry) ? lines[0] : lines[1]
    def liverLine = lines[0] == tumorLine ? lines[1] : lines[0]
    def midline = createMidlineString(liverLine.ROI.geometry, tumorLine.ROI.geometry)

    if (!midline.intersects(capsuleGeometry)) {
        addDebugAnnotations(tumorLine, liverLine, midline, capsule)
        throw new LocalRuntimeException('Expected (midline.intersects(capsuleGeometry))')
    }
    def capsuleParts = GeometryTools.splitGeometryByLineStrings(capsuleGeometry, [midline])
    if (capsuleParts.size() < 2) {
        addDebugAnnotations(tumorLine, liverLine, midline, capsule)
        throw new LocalRuntimeException('Expected (capsuleParts.size() >= 2)')
    }
    def tumorCapsule = mergeGeometries(capsuleParts.findAll { it.touches(tumor.geometry) })
    if (tumorCapsule == null) {
        addDebugAnnotations(tumorLine, liverLine, midline, capsule)
        throw new LocalRuntimeException('Expected (tumorCapsule != null)')
    }
    def liverCapsule = mergeGeometries(capsuleParts.findAll { !it.touches(tumor.geometry) })
    if (liverCapsule == null) {
        addDebugAnnotations(tumorLine, liverLine, midline, capsule)
        throw new LocalRuntimeException('Expected (liverCapsule != null)')
    }
    return [midline, liverCapsule, tumorCapsule].collect { GeometryTools.geometryToROI(it, capsule.imagePlane) }
}

def addDebugAnnotations(PathObject tumorLine, PathObject liverLine, Geometry midline, ROI capsule) {
    List<PathObject> debugAnnotations = []
    debugAnnotations << (getAnnotation(tumorLine.ROI, "00_debug_tumorLine", makeRGB(255, 50, 50)))
    debugAnnotations << (getAnnotation(liverLine.ROI, "00_debug_liverLine", makeRGB(255, 50, 50)))
    debugAnnotations << (getAnnotation(GeometryTools.geometryToROI(midline, capsule.imagePlane), "00_debug_midline", makeRGB(255, 50, 50)))
    debugAnnotations << (getAnnotation(capsule, "00_debug_capsuleGeometry", makeRGB(255, 50, 50)))
    addObjects(debugAnnotations)
}

Geometry createMidlineString(Geometry line1, Geometry line2) {
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
    def line1EndPoint = line1.getEndPoint()
    if (line2StartPoint.distance(line1StartPoint) > line2EndPoint.distance(line1StartPoint)) {
        line2StartPoint = line2EndPoint
        line2EndPoint = line2.getStartPoint()
    }
    def refLineStart = new Coordinate((line1StartPoint.x + line2StartPoint.x) / 2.0, (line1StartPoint.y + line2StartPoint.y) / 2.0)
    def refLineEnd = new Coordinate((line1EndPoint.x + line2EndPoint.x) / 2.0, (line1EndPoint.y + line2EndPoint.y) / 2.0)
    def referenceLine = geomFactory.createLineString([refLineStart, refLineEnd] as Coordinate[])
    def xDirRef = 1
    def yDirRef = (refLineStart.y - refLineEnd.y) / (refLineStart.x - refLineEnd.x)
    def xDirOrthogonal = -yDirRef
    def yDirOrthogonal = xDirRef
    def refLineLength = referenceLine.getLength()

    // Then we create a LengthIndexedLine from the reference line

    def lengthIndexedRefLine = new LengthIndexedLine(referenceLine)
    def MAX_POWER = 5
    for (int midLineResolutionPower = 0; midLineResolutionPower < MAX_POWER; midLineResolutionPower++) {
        def sampleCount = Math.max(line1.numPoints, line2.numPoints) * (2 ^ midLineResolutionPower)
        print('Trying to create a capsule midline with resolution power ' + midLineResolutionPower + ', sampleCount ' + sampleCount)
        def midpoints = [refLineStart]

        // We loop over the sampleCount,
        for (int i = 0; i <= sampleCount; i++) {
            def lengthIndex = i * (refLineLength / sampleCount)
            def referencePoint = lengthIndexedRefLine.extractPoint(lengthIndex)
            // We find the intersection between a line going through the reference point and the line2 and that is orthogonal to the reference line
            def firstOrthStart = new Coordinate(referencePoint.x - xDirOrthogonal * refLineLength / 2, referencePoint.y - yDirOrthogonal * refLineLength / 2)
            def firstOrthEnd = new Coordinate(referencePoint.x + xDirOrthogonal * refLineLength / 2, referencePoint.y + yDirOrthogonal * refLineLength / 2)
            def orthogonalLineThroughReference = geomFactory.createLineString([firstOrthStart, firstOrthEnd] as Coordinate[])
            def p1 = line1.intersection(orthogonalLineThroughReference)
            def p2 = line2.intersection(orthogonalLineThroughReference)
            if (p1 == null || p2 == null || p1.isEmpty() || p2.isEmpty()) {
                // This is the case of orthogonal line on the edges that is not intersecting one of the 2 lines, this is normal
//                print('Skipping point ' + i + ' of ' + sampleCount + ', no intersection found')
                continue
            }
            // Calculate the midpoint
            midpoints << new Coordinate((p1.centroid.x + p2.centroid.x) / 2.0, (p1.centroid.y + p2.centroid.y) / 2.0)
        }
        midpoints << refLineEnd

        def midLine = geomFactory.createLineString(midpoints as Coordinate[])
        if (!midLine.intersects(line1) && !midLine.intersects(line2)) {
            return midLine
        }
        // If the midline intersects the lines, we need to try again with a higher resolution
    }
    throw new LocalRuntimeException("Could not find a midline between the two lines that is not intersecting them")
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