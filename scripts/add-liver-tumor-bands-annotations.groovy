//file:noinspection GrMethodMayBeStatic


import groovy.transform.ImmutableOptions
import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.Geometry
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.roi.interfaces.ROI
import qupath.lib.roi.GeometryTools

import static qupath.lib.gui.scripting.QPEx.*

print('START: main')
main()
print('END: main')

def main() {
    List<Integer> liverBands = [100] * 6 + [500] * 5
    List<Integer> tumorBands = [100] * 6 + [500] * 5
    cleanupAutoAnnotations()
    List<TissueWithLines> tissuesWithLine = findTissueWithTumorLines()
    Integer coreIndex = 0
    for (TissueWithLines tissueAndLines : tissuesWithLine) {
        def (liver, tumor, capsule) = getSeparatedTissueParts(tissueAndLines)
        def liverWC = unionROI(liver, capsule)
        def tumorWC = unionROI(tumor, capsule)
        List<PathObject> tissueAnnotations = []
//        tissueAnnotations << getAnnotation(capsule, "${coreIndex}_whole_capsule", makeRGB(0, 150, 0))
//        def (liverCapsule, tumorCapsule) = splitROIInHalves(capsule, liver, tumor)
//        tissueAnnotations << getAnnotation(tumorCapsule, "${coreIndex}_tumor_capsule", makeRGB(0, 150, 0))
//        tissueAnnotations << getAnnotation(liverCapsule, "${coreIndex}_liver_capsule", makeRGB(0, 150, 0))
        tissueAnnotations << getAnnotation(liver, "${coreIndex}_liver", makeRGB(150, 150, 0))
        def (biggestLiverExpansion, liverExpansionAnnotations) = annotateHalfWithExpansions("${coreIndex}_tumor", liverWC, tumorWC, liverBands)
        tissueAnnotations.addAll(liverExpansionAnnotations)
        tissueAnnotations << getAnnotation(createCentralROI(tumor, biggestLiverExpansion), "${coreIndex}_tumor_central", makeRGB(0, 150, 0))

        tissueAnnotations << getAnnotation(tumor, "${coreIndex}_tumor", makeRGB(150, 150, 0))
        def (biggestTumorExpansion, tumorExpansionAnnotations) = annotateHalfWithExpansions("${coreIndex}_liver", tumorWC, liverWC, tumorBands)
        tissueAnnotations.addAll(tumorExpansionAnnotations)
        tissueAnnotations << getAnnotation(createCentralROI(liver, biggestTumorExpansion), "${coreIndex}_liver_central", makeRGB(0, 150, 0))
        addObjects(tissueAnnotations)
        coreIndex++
    }

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
        throw new Error("More than 2 line annotations found for tissue [${tissue?.getID()}]")
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
    def tumor = halvesGeometries.find { isTumor(it) }
    if (tumor == null) {
        throw new Error("No tumor geometry found in tissue [${tissue?.getID()}], please check the annotations")
    }
    if (halvesGeometries.size() == 2) {
        def liver = halvesGeometries.find { it != tumor }
        return [liver, tumor].collect { GeometryTools.geometryToROI(it, tissue.ROI.imagePlane) }
    }
    if (halvesGeometries.size() !== 3) {
        throw new Error("Expected 2 or 3 halves, but got ${halvesGeometries.size()}")
    }
    def geometriesTouchingTumor = halvesGeometries.findAll { it.touches(tumor) }
    if (geometriesTouchingTumor.size() != 1) {
        throw new Error("Expected 1 geometry that touches the tumor geometry, but got ${geometriesTouchingTumor.size()}")
    }
    def (capsule) = geometriesTouchingTumor

    def liver = halvesGeometries.find { it != tumor && it != capsule }
    return [liver, tumor, capsule].collect { GeometryTools.geometryToROI(it, tissue.ROI.imagePlane) }
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

Tuple<ROI> splitROIInHalves(ROI capsule, ROI liver, ROI tumor) {
    def capsuleGeometry = capsule.geometry
    def tumorGeometry = tumor.geometry
    def liverGeometry = liver.geometry
    def midline = createMidlineString(capsuleGeometry, tumorGeometry, liverGeometry)
    assert midline.intersects(capsuleGeometry)
    def (left, right) = GeometryTools.splitGeometryByLineStrings(capsuleGeometry, [midline])
    assert left != null
    assert right != null
    def tumorCapsule = left.touches(tumorGeometry) ? left : right
    def liverCapsule = left === tumorCapsule ? right : left
    return [liverCapsule, tumorCapsule].collect { GeometryTools.geometryToROI(it, capsule.imagePlane) }
}

Geometry createMidlineString(Geometry geometry, Geometry tumorGeometry, Geometry liverGeometry) {
    def liverFacingWall = geometry.intersection(liverGeometry)
    def tumorFacingWall = geometry.intersection(tumorGeometry)
    return createMidlineIn(liverFacingWall, tumorFacingWall, geometry)
}

Geometry createMidlineIn(Geometry line1, Geometry line2, Geometry shapeThatShouldBeSplit) {
    def factory = shapeThatShouldBeSplit.getFactory()
    def boundary = shapeThatShouldBeSplit.getBoundary()

    // Get boundary intersection points
    def inter1 = line1.intersection(boundary)
    def inter2 = line2.intersection(boundary)

    if (inter1.isEmpty() || inter2.isEmpty()) {
        throw new Error("Failed to find boundary intersections.")
    }

    def A = inter1.getCoordinates()[0]
    def B = inter2.getCoordinates()[0]

    // Convert lines to coordinate arrays
    def coords1 = line1.getCoordinates()
    def coords2 = line2.getCoordinates()

    // Find starting points that are closest to each other
    def start1 = coords1.min { it.distance(A) }
    def start2 = coords2.min { it.distance(B) }

    // Walk along both lines simultaneously to find corresponding points
    def midPoints = []
    def idx1 = coords1.findIndexOf {it === start1}
    def idx2 = coords2.findIndexOf {it === start2}

    // Determine walking directions
    def dir1 = A.distance(coords1[0]) < A.distance(coords1[-1]) ? 1 : -1
    def dir2 = B.distance(coords2[0]) < B.distance(coords2[-1]) ? 1 : -1

    // Walk until we reach the end of either line
    while (idx1 >= 0 && idx1 < coords1.size() && idx2 >= 0 && idx2 < coords2.size()) {
        def p1 = coords1[idx1]
        def p2 = coords2[idx2]

        // Add midpoint
        midPoints << new Coordinate((p1.x + p2.x)/2, (p1.y + p2.y)/2)

        // Move forward along both lines
        idx1 += dir1
        idx2 += dir2

        // Optional: Skip some points for efficiency if lines are very dense
        if (coords1.size() > 100 || coords2.size() > 100) {
            idx1 += (coords1.size()/50).toInteger()
            idx2 += (coords2.size()/50).toInteger()
        }
    }

    // Sort points by progression along the path
    def allPoints = [A] + sortPointsAlongPath(midPoints) + [B]

    // Create smooth path
    def smoothedPath = smoothPath(allPoints, 3)

    return factory.createLineString(smoothedPath as Coordinate[])
}

List<Coordinate> sortPointsAlongPath(List<Coordinate> points) {
    if (points.size() <= 2) return points

    def sorted = [points[0]]
    def remaining = points[1..-1] as List

    while (!remaining.isEmpty()) {
        def last = sorted.last()
        def closest = remaining.min { it.distance(last) }
        sorted << closest
        remaining.remove(closest)
    }

    return sorted
}

List<Coordinate> smoothPath(List<Coordinate> path, int iterations = 1) {
    if (path.size() <= 2) return path

    def smoothed = path.clone()

    iterations.times {
        def newPath = [smoothed[0]]
        for (int i = 1; i < smoothed.size()-1; i++) {
            // Stronger smoothing towards neighbors
            newPath << new Coordinate(
                    (smoothed[i-1].x * 0.4 + smoothed[i].x * 0.2 + smoothed[i+1].x * 0.4),
                    (smoothed[i-1].y * 0.4 + smoothed[i].y * 0.2 + smoothed[i+1].y * 0.4)
            )
        }
        newPath << smoothed[-1]
        smoothed = newPath
    }

    return smoothed
}