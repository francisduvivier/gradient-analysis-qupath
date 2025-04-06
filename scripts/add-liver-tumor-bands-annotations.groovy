import org.locationtech.jts.geom.Geometry
import qupath.lib.common.GeneralTools
import qupath.lib.geom.Point2
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.roi.GeometryTools
import qupath.lib.roi.ROIs
import static qupath.lib.gui.scripting.QPEx.* // For getAnnotationObjects(), addObject()

def annotations = getAnnotationObjects()

// Find required annotations
def tissue = annotations.find { it.getROI().isArea() }
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

    def tissuePoints = tissue.getROI().getAllPoints().subList(0, 1000)
    def polygonROI = ROIs.createPolygonROI(tumorPoints+tissuePoints, plane)

    // Verify the created ROI is a valid area
    if (!polygonROI.isArea()) {
        print "Created shape is not a valid area - check tumor line configuration"
        return
    }

    // Create and a
    def tumorArea = PathObjects.createAnnotationObject(polygonROI)
    tumorArea.setName("Tumor Boundary")
    addObject(tumorArea)
    print "Successfully added tumor boundary annotation"

} catch (Exception e) {
    print "Could not create tumor area: " + e.getMessage()
}