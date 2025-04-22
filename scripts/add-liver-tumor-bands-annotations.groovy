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
        def (liverCapsule, tumorCapsule) = splitROIInHalves(capsule, liver, tumor)
//        tissueAnnotations << getAnnotation(capsule, "${coreIndex}_whole_capsule", makeRGB(0, 150, 0))
        tissueAnnotations << getAnnotation(tumorCapsule, "${coreIndex}_tumor_capsule", makeRGB(0, 150, 0))
        tissueAnnotations << getAnnotation(liverCapsule, "${coreIndex}_liver_capsule", makeRGB(0, 150, 0))
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

    def inter1 = line1.intersection(boundary)
    def inter2 = line2.intersection(boundary)

    if (inter1.isEmpty() || inter2.isEmpty()) {
        throw new Error("Failed to find boundary intersections.")
    }

    def A = inter1.getCoordinates()[0]
    def B = inter2.getCoordinates()[0]

    // Compute the midpoint between the boundary intersections
    def mid = new Coordinate((A.x + B.x) / 2.0, (A.y + B.y) / 2.0)

    return factory.createLineString([A, mid, B] as Coordinate[])
}