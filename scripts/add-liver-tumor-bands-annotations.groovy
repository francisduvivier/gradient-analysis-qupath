import org.locationtech.jts.geom.Geometry
import qupath.lib.geom.Point2
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.GeometryROI
import qupath.lib.roi.PolygonROI
import qupath.lib.roi.ROIs

import static qupath.lib.gui.scripting.QPEx.*

// For getAnnotationObjects(), addObject()

def annotations = getAnnotationObjects()

// Find required annotations
def tissue = annotations.find { it.getROI().isArea() && it.getName() == 'mTissue' }
def tumorLine = annotations.find { it.getROI().isLine() }

if (tissue == null || tumorLine == null) {
    String missing = []
    if (tissue == null) missing.add("a closed tissue annotation")
    if (tumorLine == null) missing.add("a tumor line annotation")
    print "Missing required annotations: need " + missing.join(" and ") + "."
    return
}

print "Found tissue annotation: ${tissue.getName()}"
print "Found tumor line annotation: ${tumorLine.getName()}"

// Get the tumor line points - corrected method call
def tumorLineROI = tumorLine.getROI()
def tumorPoints = tumorLineROI.getAllPoints()  // Correct method name for PolylineROI

// Validate we can create a polygon
if (tumorPoints.size() < 3) {
    print "Cannot create area - tumor line needs at least 3 points (has ${tumorPoints.size()})"
    return
}

try {
    def halves = getSeparatedTissuePoints(tissue, tumorLine)
    addHalfTissueAnnotation(halves[0], 'half[0]')
    addHalfTissueAnnotation(halves[1], 'half[1]')

} catch (Exception e) {
    print "Could not create tumor area: " + e.getMessage()
}

private void addHalfTissueAnnotation(PolygonROI halfTissueROI, String name) {
    def halfTissue = PathObjects.createAnnotationObject(halfTissueROI)
    halfTissue.setName(name)
    addObject(halfTissue)
    print "Successfully added halfTissue [$name] boundary annotation"
}

Tuple<PolygonROI> getSeparatedTissuePoints(PathObject tissue, PathObject roughTumorLine) {
    def plane = tissue.getROI().getImagePlane()

    Geometry interSection = roughTumorLine.getROI().getGeometry().intersection(tissue.getROI().getGeometry())
    def tissuePoints = tissue.getROI().getAllPoints()
    def tumorLinePoints = interSection.getCoordinates().collect { new Point2(it.x, it.y) }

    def closestToStartOfTumorIndex = findClosestPointIndex(tumorLinePoints.first, tissuePoints)
    print "closestToStartOfTumorIndex: " + closestToStartOfTumorIndex
    def closestToEndOfTumorIndex = findClosestPointIndex(tumorLinePoints.last, tissuePoints)
    print "closestToEndOfTumorIndex: " + closestToEndOfTumorIndex
    def closestPointStartIndex = (int) Math.min(closestToStartOfTumorIndex, closestToEndOfTumorIndex)
    def closestPointEndIndex = (int) Math.max(closestToStartOfTumorIndex, closestToEndOfTumorIndex)
    List<Point2> halfInInnerPartOfArray = tissuePoints.subList(closestPointStartIndex, closestPointEndIndex)
    print "topTissuePoints: " + halfInInnerPartOfArray.size()

    def halfInOuterPartOfArray = tissuePoints.subList(closestPointEndIndex, tissuePoints.size()) + tissuePoints.subList(0, closestPointStartIndex)

    return [halfInInnerPartOfArray, halfInOuterPartOfArray].collect { combineLinesToRoi(tumorLinePoints, it, plane) }
}

static int findClosestPointIndex(Point2 point, List<Point2> otherPoints) {
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

static PolygonROI combineLinesToRoi(List<Point2> firstPoints, List<Point2> lastPoints, ImagePlane plane) {
    if (firstPoints.last.distance(lastPoints.first) > firstPoints.last.distance(lastPoints.last)) {
        // If the last point of the separator does not closely connect to the first point of the tissue, we need to connect them the other way around
        lastPoints = lastPoints.reversed()
    }
    return ROIs.createPolygonROI(firstPoints + lastPoints, plane)
}