//file:noinspection GrMethodMayBeStatic
import org.locationtech.jts.geom.Geometry
import qupath.lib.geom.Point2
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.PolygonROI
import qupath.lib.roi.ROIs
import qupath.lib.roi.interfaces.ROI

import static qupath.lib.gui.scripting.QPEx.*

print('START: main')
main()
print('END: main')

def main() {
    List<Tuple<PathObject>> tissuesWithLine = findTissueWithTumorLines()
    Integer coreIndex = 0
    for (Tuple<PathObject> tissueAndLine : tissuesWithLine) {
        def (tissue, tumorLine) = tissueAndLine
        def halves = getSeparatedTissuePoints(tissue, tumorLine)
        annotateHalfWithExpansions("liver[$coreIndex]", halves[1], halves[0], 4, 100)
        annotateHalfWithExpansions("tumor[$coreIndex]", halves[0], halves[1], 6, 100)
        coreIndex++
    }

}


void annotateHalfWithExpansions(String name, PolygonROI polygonROI, PolygonROI otherHalf, int amount, int microns) {
    addAnnotation(polygonROI, name, makeRGB(150, 150, 0))
    double distance = getDistance(microns)
    def prevExpansion = polygonROI.geometry
    for (int i = 1; i <= amount; i++) {
        def expansion = polygonROI.geometry.buffer(i * distance)
        def expansionBand = expansion.difference(prevExpansion)
        prevExpansion = expansion
        def intersection = expansionBand.intersection(otherHalf.geometry)

        def bandColor = makeRGB(20 * i, 40 * i, 200 - 30 * i)
        addAnnotation(getROIForGeometry(intersection, polygonROI.imagePlane), "$name [${i * microns} micrometer]", bandColor)
    }
}

List<Tuple<PathObject>> findTissueWithTumorLines() {
    Collection<PathObject> annotations = getAnnotationObjects()
    // Find required annotations
    return annotations.collect { getTissueWithLine(it) }.findAll { it !== null }
}

Tuple<PathObject> getTissueWithLine(PathObject candidate) {
    if (!candidate.ROI.isArea()) {
        return null
    }

    def tissue = candidate
    print "Found tissue annotation: ${tissue?.name}"

    Collection<PathObject> annotations = getAnnotationObjects()
    def tumorLine = annotations.find { it.ROI.isLine() && it.ROI.geometry.intersects(tissue.ROI.geometry) }
    print "Found tumor line annotation: ${tumorLine?.name}"

    if (tissue == null || tumorLine == null) {
        String missing = []
        if (tissue == null) missing.add("a closed tissue annotation")
        if (tumorLine == null) missing.add("a tumor line annotation")
        throw new RuntimeException("Missing required annotations: need " + missing.join(" and ") + ".")
    }
    return [tissue, tumorLine]
}

PathObject addAnnotation(ROI roi, String name, Integer color = makeRGB(255, 255, 0)) {
    def halfTissue = PathObjects.createAnnotationObject(roi)
    halfTissue.setColor(color)
    halfTissue.setName(name)
    Collection<PathObject> annotations = getAnnotationObjects()
    def existing = annotations.find { it.name == name }
    if (existing != null) {
        print "Removing existing annotation [$name]"
        removeObject(existing, false)
    }
    addObject(halfTissue)
    print "Successfully added ROI [$name] annotation"
    return halfTissue
}

Tuple<PolygonROI> getSeparatedTissuePoints(PathObject tissue, PathObject roughTumorLine) {
    def plane = tissue.ROI.imagePlane

    Geometry interSection = roughTumorLine.ROI.getGeometry().intersection(tissue.ROI.getGeometry())
    def tissuePoints = tissue.ROI.getAllPoints()
    def tumorLinePoints = getGeometryPoints(interSection)

    def closestToStartOfTumorIndex = findClosestPointIndex(tumorLinePoints.first, tissuePoints)
    print "closestToStartOfTumorIndex: " + closestToStartOfTumorIndex
    def closestToEndOfTumorIndex = findClosestPointIndex(tumorLinePoints.last, tissuePoints)
    print "closestToEndOfTumorIndex: " + closestToEndOfTumorIndex
    def closestPointStartIndex = (int) Math.min(closestToStartOfTumorIndex, closestToEndOfTumorIndex)
    def closestPointEndIndex = (int) Math.max(closestToStartOfTumorIndex, closestToEndOfTumorIndex)
    List<Point2> halfInInnerPartOfArray = tissuePoints.subList(closestPointStartIndex, closestPointEndIndex + 1)
    print "topTissuePoints: " + halfInInnerPartOfArray.size()

    def halfInOuterPartOfArray = tissuePoints.subList(closestPointEndIndex, tissuePoints.size()) + tissuePoints.subList(0, closestPointStartIndex + 1)

    return [halfInInnerPartOfArray, halfInOuterPartOfArray].collect { combineLinesToRoi(tumorLinePoints, it, plane) }
}


int findClosestPointIndex(Point2 point, List<Point2> otherPoints) {
    def closestIndex = 0
    def closestDistance = otherPoints[closestIndex].distance(point)
    for (int i = 1; i < otherPoints.size(); i++) {
        def otherPoint = otherPoints[i]
        def distance = otherPoint.distance(point)
        if (distance < closestDistance) {
            closestIndex = i
            closestDistance = distance
        }
    }
    return closestIndex
}

PolygonROI combineLinesToRoi(List<Point2> firstPoints, List<Point2> lastPoints, ImagePlane plane) {
    if (firstPoints.last.distance(lastPoints.first) > firstPoints.last.distance(lastPoints.last)) {
        // If the last point of the separator does not closely connect to the first point of the tissue, we need to connect them the other way around
        lastPoints = lastPoints.reversed()
    }
    return ROIs.createPolygonROI(firstPoints + lastPoints, plane)
}


ROI getROIForGeometry(Geometry geometry, ImagePlane plane) {
    return ROIs.createPolygonROI(getGeometryPoints(geometry), plane)
}

List<Point2> getGeometryPoints(Geometry interSection) {
    interSection.getCoordinates().collect { new Point2(it.x, it.y) }
}

double getDistance(double microns) {
    def imageData = getCurrentImageData()
    def server = imageData.getServer()
// We need the pixel size
    def cal = server.getPixelCalibration()
    if (!cal.hasPixelSizeMicrons()) {
        print 'We need the pixel size information here!'
        throw new RuntimeException("We need the pixel size information here!")
    }
    return microns / cal.getAveragedPixelSizeMicrons()
}