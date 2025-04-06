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
    def (PathObject tissue, tumorLine) = findTissueWithTumorLine()
    def halves = getSeparatedTissuePoints(tissue, tumorLine)
    addAnnotation(halves[0], 'half[0]')
    addAnnotation(halves[1], 'half[1]')
}

Tuple<PathObject> findTissueWithTumorLine() {
    Collection<PathObject> annotations = getAnnotationObjects()

    // Find required annotations
    def tissue = annotations.find { it.ROI.isArea() && it.getName() == 'mTissue' }
    print "Found tissue annotation: ${tissue?.name}"
    def tumorLine = annotations.find { it.ROI.isLine() && it.ROI.geometry.intersects(tissue.ROI.geometry) }
    print "Found tumor line annotation: ${tumorLine?.name}"

    if (tissue == null || tumorLine == null) {
        String missing = []
        if (tissue == null) missing.add("a closed tissue annotation")
        if (tumorLine == null) missing.add("a tumor line annotation")
        throw new RuntimeException("Missing required annotations: need " + missing.join(" and ") + ".")
    }
    [tissue, tumorLine]
}

void addAnnotation(ROI roi, String name) {
    def halfTissue = PathObjects.createAnnotationObject(roi)
    halfTissue.setName(name)
    addObject(halfTissue)
    print "Successfully added ROI [$name] annotation"
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

List<Point2> getGeometryPoints(Geometry interSection) {
    interSection.getCoordinates().collect { new Point2(it.x, it.y) }
}
