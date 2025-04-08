//file:noinspection GrMethodMayBeStatic
import org.locationtech.jts.geom.Geometry
import qupath.lib.geom.Point2
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.PolygonROI
import qupath.lib.roi.ROIs
import qupath.lib.roi.interfaces.ROI
import qupath.lib.roi.GeometryTools

import static qupath.lib.gui.scripting.QPEx.*

print('START: main')
main()
print('END: main')

def main() {
//    cleanupAutoAnnotations()
    List<Tuple<PathObject>> tissuesWithLine = findTissueWithTumorLines()
    Integer coreIndex = 0
    for (Tuple<PathObject> tissueAndLine : tissuesWithLine) {
        def (tissue, tumorLine) = tissueAndLine
        def (liver, tumor) = getSeparatedTissuePoints(tissue, tumorLine)

        addAnnotation(liver, "${coreIndex}_liver", makeRGB(150, 150, 0))
        def biggestLiverExpansion = annotateHalfWithExpansions("${coreIndex}_tumor", liver, tumor, 6, 100)
        addAnnotation(createCentralROI(tumor, biggestLiverExpansion), "${coreIndex}_tumor_central", makeRGB(0, 150, 0))

        addAnnotation(tumor, "${coreIndex}_tumor", makeRGB(150, 150, 0))
        def biggestTumorExpansion = annotateHalfWithExpansions("${coreIndex}_liver", tumor, liver, 4, 100)
        addAnnotation(createCentralROI(liver, biggestTumorExpansion), "${coreIndex}_liver_central", makeRGB(0, 150, 0))

        coreIndex++
    }

}


Geometry annotateHalfWithExpansions(String name, ROI startPolygon, ROI intersectingPolygon, int amount, int microns) {
    double distance = getDistance(microns)
    def prevExpansion = startPolygon.geometry
    for (int i = 1; i <= amount; i++) {
        def expansion = startPolygon.geometry.buffer(i * distance)
        def expansionBand = expansion.difference(prevExpansion)
        prevExpansion = expansion
        def intersection = expansionBand.intersection(intersectingPolygon.geometry)

        def bandColor = makeRGB(20 * i, 40 * i, 200 - 30 * i)
        addAnnotation(GeometryTools.geometryToROI(intersection, startPolygon.imagePlane), "${name}_${i * microns}µm", bandColor)
    }
    return prevExpansion
}

List<Tuple<PathObject>> findTissueWithTumorLines() {
    Collection<PathObject> annotations = getAnnotationObjects()
    // Find required annotations
    return annotations.collect { getTissueWithLine(it) }.findAll { it !== null }
}

Tuple<PathObject> getTissueWithLine(PathObject candidate) {
    if (!candidate.ROI.isArea() || candidate.classifications.contains('Auto')) {
        return null
    }

    def tissue = candidate
    print "Found tissue annotation: ${tissue?.name}"

    Collection<PathObject> annotations = getAnnotationObjects()
    def tumorLine = annotations.find { it.ROI.isLine() && it.ROI.geometry.intersects(tissue.ROI.geometry) }
    print "Found tumor line annotation: ${tumorLine?.name}"

    if (tissue == null || tumorLine == null) {
        String[] missing = []
        if (tissue == null) missing.add("a closed tissue annotation")
        if (tumorLine == null) missing.add("a tumor line annotation")
        throw new RuntimeException("Missing required annotations: need " + missing.join(" and ") + ".")
    }
    return [tissue, tumorLine]
}

PathObject addAnnotation(ROI roi, String name, Integer color = makeRGB(255, 255, 0)) {
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
    addObject(newAnnotation)
    print "Successfully added ROI [$name] annotation"
    return newAnnotation
}

Tuple<ROI> getSeparatedTissuePoints(PathObject tissue, PathObject roughTumorLine) {
    def halvesGeometries = GeometryTools.splitGeometryByLineStrings(tissue.ROI.geometry, [roughTumorLine.ROI.geometry])
    def halves = halvesGeometries.collect { GeometryTools.geometryToROI(it, tissue.ROI.imagePlane)}
    if (isTumor(halves[0])) {
        halves = halves.reverse()
    }
    return halves
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

boolean isTumor(ROI polygonROI) {
    Collection<PathObject> annotations = getAnnotationObjects()
    return annotations.find { it.classifications.contains('Tumor') && polygonROI.geometry.intersects(it.ROI.geometry)} != null
}

ROI createCentralROI(ROI startRoi, Geometry biggestExpansion) {
    def centralGeometry = startRoi.geometry.difference(biggestExpansion)
    return GeometryTools.geometryToROI(centralGeometry, startRoi.imagePlane)
}

void cleanupAutoAnnotations() {
    Collection<PathObject> annotations = getAnnotationObjects()
    for(PathObject annotation : annotations) {
        if (annotation.classifications.contains('Auto')) {
            removeObject(annotation, false)
        }
    }
}