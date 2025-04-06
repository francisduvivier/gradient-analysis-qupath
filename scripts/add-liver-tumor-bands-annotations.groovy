import org.locationtech.jts.geom.Geometry
import qupath.lib.geom.Point2
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.GeometryROI
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
    // Create polygon ROI (automatically connects last to first point)
    def plane = tumorLineROI.getImagePlane()

    Geometry interSection = tumorLineROI.getGeometry().intersection(tissue.getROI().getGeometry())
    def tumorLineWithinTissueROI = new GeometryROI(interSection, plane)

    def halves = getSeparatedTissuePoints(tissue, tumorLineWithinTissueROI)
    addHalfTissueAnnotation(tumorLineWithinTissueROI, halves[0], plane, 'half[0]')
    addHalfTissueAnnotation(tumorLineWithinTissueROI, halves[1], plane, 'half[1]')

} catch (Exception e) {
    print "Could not create tumor area: " + e.getMessage()
}

private void addHalfTissueAnnotation(GeometryROI tumorLineWithinTissueROI, List<Point2> halfTissuePoints, ImagePlane plane, String name) {
    def separatorPoints = tumorLineWithinTissueROI.getAllPoints()
    if (separatorPoints.last.distance(halfTissuePoints.first) > separatorPoints.last.distance(halfTissuePoints.last)) {
        // If the last point of the separator does not closely connect to the first point of the tissue, we need to connect them the other way around
        halfTissuePoints = halfTissuePoints.reversed()
    }

    def topTissueROI = ROIs.createPolygonROI(separatorPoints + halfTissuePoints, plane)

    def tumorArea = PathObjects.createAnnotationObject(topTissueROI)
    tumorArea.setName(name)
    addObject(tumorArea)
    print "Successfully added [$name] boundary annotation"
}

def getSeparatedTissuePoints(PathObject tissue, GeometryROI tumorLineWithinTissueROI) {
    def tissuePoints = tissue.getROI().getAllPoints()
    def tumorLinePoints = tumorLineWithinTissueROI.getAllPoints()
    def closestToStartOfTumorIndex = findClosestPointIndex(tumorLinePoints.first, tissuePoints)
    print "closestToStartOfTumorIndex: " + closestToStartOfTumorIndex
    def closestToEndOfTumorIndex = findClosestPointIndex(tumorLinePoints.last, tissuePoints)
    print "closestToEndOfTumorIndex: " + closestToEndOfTumorIndex
    def closestPointStartIndex = (int) Math.min(closestToStartOfTumorIndex, closestToEndOfTumorIndex)
    def closestPointEndIndex = (int) Math.max(closestToStartOfTumorIndex, closestToEndOfTumorIndex)
    List<Point2> topTissuePoints = tissuePoints.subList(closestPointStartIndex, closestPointEndIndex)
    print "topTissuePoints: " + topTissuePoints.size()

    def halfInOuterPartOfArray = tissuePoints.subList(closestPointEndIndex, tissuePoints.size())+ tissuePoints.subList(0, closestPointStartIndex)
    return [topTissuePoints, halfInOuterPartOfArray]
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