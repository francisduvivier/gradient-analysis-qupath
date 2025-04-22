//file:noinspection GrMethodMayBeStatic


import groovy.transform.ImmutableOptions
import groovy.transform.ImmutableProperties
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
        def (liverWC, tumorWC) = getSeparatedTissueParts(tissueAndLines)
        def tumor = diffRoi(tumorWC, liverWC)
        def liver = diffRoi(liverWC, tumorWC)
        List<PathObject> tissueAnnotations = []
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
record TissueWithLines(PathObject tissue, Collection<PathObject> lines){}

TissueWithLines getTissueWithLine(PathObject candidate) {
    assert candidate != null
    if (!candidate.ROI.isArea() || candidate.classifications.contains('Auto')) {
        return null
    }

    def tissue = candidate
    print "Found tissue annotation: ${tissue?.name}"

    Collection<PathObject> annotations = getAnnotationObjects()
    Collection<PathObject> tissueLines = findLinesForTissue(annotations, tissue)

    return new TissueWithLines(tissue, tissueLines)
}

Collection<PathObject> findLinesForTissue(Collection<PathObject> annotations, PathObject tissue) {
    Collection<PathObject> lines = annotations.findAll { it.ROI.isLine() && it.ROI.geometry.intersects(tissue.ROI.geometry) }
    if (lines.size() == 0) {
        throw new Error("No tumor line annotation found for tissue [${tissue?.getID()}]")
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
    def halvesGeometries = GeometryTools.splitGeometryByLineStrings(tissue.ROI.geometry, tumorLines.collect {it.ROI.geometry})
    if(halvesGeometries.size() == 3) {
        print("found tree parts for tissue [${tissue.getID()}], merging to make overlapping parts")
        halvesGeometries = [halvesGeometries[0].union(halvesGeometries[2]), halvesGeometries[1].union(halvesGeometries[2])]
    }
    if(halvesGeometries.size() > 3) {
        throw new Error("Expected 2 or 3 halves, but got ${halvesGeometries.size()}")
    }
    List<ROI> halves = halvesGeometries.collect { GeometryTools.geometryToROI(it, tissue.ROI.imagePlane) }
    if (isTumor(halves[0])) {
        halves = halves.reverse()
    }
    def (liverWC, tumorWC) = halves
    return [liverWC, tumorWC]
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
    return annotations.find { it.classifications.contains('Tumor') && polygonROI.geometry.intersects(it.ROI.geometry) } != null
}

ROI createCentralROI(ROI startRoi, Geometry biggestExpansion) {
    def centralGeometry = startRoi.geometry.difference(biggestExpansion)
    return GeometryTools.geometryToROI(centralGeometry, startRoi.imagePlane)
}

void cleanupAutoAnnotations() {
    Collection<PathObject> annotations = getAnnotationObjects()
    removeObjects(annotations.findAll { it.classifications.contains('Auto') }, false)
}

ROI diffRoi(ROI roi1, ROI roi2) {
    return GeometryTools.geometryToROI(roi1.geometry.difference(roi2.geometry), roi1.imagePlane)
}