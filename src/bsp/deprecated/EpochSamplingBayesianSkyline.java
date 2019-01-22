package bsp.deprecated;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;
import bsp.distributions.SegmentType;

import java.util.*;

/**
 * Epoch Sampling Bayesian Skyline plot
 *
 * Population size estimates are as in the regular Bayesian Skyline plot
 *
 * Linear Sampling Intensity is piecewise constant on a set of epochs
 *
 * Epoch times
 *  - Epoch times are always specified as height from the tips of the tree (age from the most-recent sample)
 *  - epoch.start < epoch.end (end is always closer to the tMRCA)
 *  - Epoch is defined on (epoch.start, epoch.end] except for the first epoch which is [0, epoch.end]
 *
 *
 */
public class EpochSamplingBayesianSkyline extends TreeDistribution {

    final public Input<RealParameter> popSizeInput =
            new Input<>("popSizes","Effective population size (skyline parameter)", Input.Validate.REQUIRED);

    final public Input<IntegerParameter> groupSizeInput =
            new Input<>("groupSizes", "The number of events in each group in the skyline (use robust design if not provided)");

    final public Input<RealParameter> samplingIntensityInput =
            new Input<>("samplingIntensity","Constant sampling intensity", Input.Validate.REQUIRED);

    final public Input<RealParameter> samplingIntensityChangeTimesInput =
            new Input<>("samplingEpochTimes","Times when the sampling epochs change (distance from most recent tip)");

    final public Input<Double> minGroupWidthInput =
            new Input<>("minGroupWidth","Minimum width for a group (end-start)",0.0);


    protected TreeInterface tree;
    protected RealParameter popSizes, samplingIntensity, epochTimes;
    protected IntegerParameter groupSizes;
    protected Double minGroupWidth;

    protected int [] cumulativeGroupSizes,
                     groupLastSegments,
                     epochLastSegments,
                     epochLastSamples,
                     epochFirstSamples;

//boolean arraysUpdated = false;
    boolean segmentsUpdated = false;

    // Segments of the skyline model
    List<EpochSamplingSkylineSegment> segments;


    @Override
    public void initAndValidate() {

        int nrEvents, nrGroups, nrEpochs;

        //////////////////////////
        // Get skyline parameter
        popSizes = popSizeInput.get();
        nrGroups = popSizes.getDimension();

        /////////////////////////////////////
        // Get sampling intensity and epochs
        samplingIntensity = samplingIntensityInput.get();
        if (samplingIntensityChangeTimesInput.get() == null) {
            epochTimes = new RealParameter(new Double [] {});
        } else {
            epochTimes = samplingIntensityChangeTimesInput.get();
        }

        nrEpochs = epochTimes.getDimension();

        if (samplingIntensity.getDimension() == 1 && nrEpochs > 1) {
            samplingIntensity.setDimension(nrEpochs);
            for (int i = 1; i < nrEpochs; i++)
                samplingIntensity.setValue(i, samplingIntensity.getValue(0));
        } else
        if (nrEpochs != samplingIntensity.getDimension()-1) {
            throw new IllegalArgumentException("Dimension of sampling intensity must be either 1 or equal to the number of sampling epochs minus 1!");
        }


        ///////////////////////
        // Get tree
        if (treeInput.get() == null) {
            throw new IllegalArgumentException("Only tree input (not tree intervals) should be specified!");
        } else {
            tree = treeInput.get();
        }
        nrEvents = tree.getNodeCount();

        ////////////////////////////
        // Get minimum group width
        minGroupWidth = minGroupWidthInput.get();

        ////////////////////
        // Get group sizes
        // If groupSizes are not specified use robust design (equal group sizes)
        if (groupSizeInput.get() != null) {
            groupSizes = groupSizeInput.get();
        } else {
            groupSizes = getRobustGroupSizes(nrEvents, nrGroups);
        }

        // Group sizes does not equal the dimension of the skyline parameter
        if (groupSizes.getDimension() != nrGroups) {
            throw new IllegalArgumentException("Number of groups should match the dimension of the effective population size.");
        }

        // More groups than coalescent events
        if (groupSizes.getDimension() > nrEvents) {
            throw new IllegalArgumentException("There are more groups than coalescent+sampling events in the tree.");
        }

        // Initalise arrays
        cumulativeGroupSizes = new int[nrGroups];
        groupLastSegments = new int [nrGroups];
        epochLastSegments = new int [nrEpochs+1];
        epochLastSamples  = new int [nrEpochs+1];
        epochFirstSamples = new int [nrEpochs+1];

        // GroupSizes needs to add up to coalescent events
        cumulativeGroupSizes[0] = groupSizes.getValue(0);
        for (int i = 1; i < cumulativeGroupSizes.length; i++) {
            cumulativeGroupSizes[i] = cumulativeGroupSizes[i-1] + groupSizes.getValue(i);
        }

        if (cumulativeGroupSizes[nrGroups - 1] != nrEvents) {
            Log.warning.println("WARNING: The sum of the initial group sizes does not match the number of coalescent " +
                                "and sampling events in the tree. Initializing to equal group sizes (robust design)");

            groupSizes.assignFromWithoutID(getRobustGroupSizes(nrEvents, nrGroups));

            // Recalculate cumulative group sizes, because group sizes have been changed
            //updateArrays();
        }


    }


    /**
     * Distribute an equal number of coalescent events to each group
     *
     * @param events
     * @param groups
     */
    protected IntegerParameter getRobustGroupSizes(int events, int groups) {

        int eventsEach   = events / groups;
        int eventsExtras = events % groups;

        Integer[] values = new Integer[groups];
        for (int i = 0; i < groups; i++) {
            if (i < eventsExtras) {
                values[i] = eventsEach + 1;
            } else {
                values[i] = eventsEach;
            }
        }

        IntegerParameter parameter = new IntegerParameter(values);
        parameter.setBounds(1, Integer.MAX_VALUE);

        return(parameter);
    }

    /**
     *
     *
     *
     */


    /**
     * Updates the arrays used in likelihood calculation and other methods
     *
     * - cumulateGroupSizes: cumulative sum of group sizes
     * - groupLastSegments:  Segment indices of the end of each group
     * - epochLastSegments:  Segment indices of the end of each epoch
     * - epochFirstSamples:  Segment indices of the first sample in each epoch
     * - epochLastSamples:   Segment indices of the last sample in each epoch
     *
     */
    protected void updateArrays() {

        // Get cumulative group sizes
        cumulativeGroupSizes[0] = groupSizes.getValue(0);
        for (int i = 1; i < cumulativeGroupSizes.length; i++) {
            cumulativeGroupSizes[i] = cumulativeGroupSizes[i-1] + groupSizes.getValue(i);
        }

        // Fill accounting arrays
        // - Segment indices of the end of each epoch
        // - Segment indices for the first and last sampling events in each epoch
        //
        // TODO: Fix code duplication
        int i = 0,              // Counts epochs
            j = 0,              // Counts segments in the tree (used for groups) - current segment is j+i
            k = 0,              // Counts groups
            nsamples    = 0,    // Counts the number of samples in an epoch
            prevSample  = -1;   // Keeps track of the segment index of the last seen sampling event
        for (EpochSamplingSkylineSegment segment : segments) {

            switch (segment.getSegmentType()) {
                case SAMPLE :
                    prevSample = j+i;

                    nsamples++;
                    if (epochFirstSamples[i] < 0)
                        epochFirstSamples[i] = j+i;

                    if (j >= cumulativeGroupSizes[k]) {
                        groupLastSegments[k] = j+i-1;
                        k++;
                    }
                    j++;

                    break;

                case COALESCENT:
                    if (j >= cumulativeGroupSizes[k]) {
                        groupLastSegments[k] = j+i-1;
                        k++;
                    }
                    j++;

                    break;

                case NOTHING:
                    //if (nsamples < 2) {
                    //    Log.warning.println("WARNING: Epoch "+i+" contains less than 2 sampling events. Model is not identifiable.");
                    //}

                    epochLastSegments[i] = j+i;
                    epochLastSamples[i]  = prevSample;

                    i++;
                    epochFirstSamples[i] = -1;
                    nsamples = 0;

                    break;

                default:
                    break;
            }

        }
        groupLastSegments[k] = j+i-1;
        epochLastSegments[i] = j+i-1;
        epochLastSamples[i]  = prevSample;


//        // Throw warnings (if some epochs have no sampling events)
//        for (int k = 1; k < epochLastSegments.length; k++) {
//            if (epochLastSegments[k-1] > epochLastSamples[k]) {
//                Log.warning.println("WARNING: Epoch "+(k-1)+" does not contain any sampling events, sampling intensity " +
//                                    "will be 0 throughout the epoch.");
//            }
//        }

//        System.out.println(Arrays.toString(cumulativeGroupSizes));
//        System.out.println(Arrays.toString(groupLastSegments));
//        System.out.println(Arrays.toString(epochLastSegments));
//        System.out.println(Arrays.toString(epochFirstSamples));
//        System.out.println(Arrays.toString(epochLastSamples));


        //arraysUpdated = true;
    }



    /**
     * 1.) Collect the end times and types of all segments
     * 2.) Orders times and types
     * 3.) Count lineages and widths of each segment
     *
     */
    protected List collectSegmentTimes() {

        Node[] nodes = tree.getNodesAsArray();
        ArrayList<EpochSamplingSkylineSegment> segments = new ArrayList<EpochSamplingSkylineSegment>();

        // Get coalescent and sampling events on the tree, epoch times and sort by end time, then segment type
        // such that for multiple segments ending at the same time, the epoch change time is the last segment
        for (int i = 0; i < nodes.length; i++) {
            segments.add(new EpochSamplingSkylineSegment(nodes[i].getHeight(),
                                                         nodes[i].isLeaf() ? SegmentType.SAMPLE : SegmentType.COALESCENT));
        }
        for (int i = 0; i < epochTimes.getDimension(); i++) {
            segments.add(new EpochSamplingSkylineSegment(epochTimes.getValue(i), SegmentType.NOTHING));
        }
        //Collections.sort(segments, (o1, o2) -> (Double.compare(o1.getEnd(),o2.getEnd())));
        segments.sort(Comparator.comparing(EpochSamplingSkylineSegment::getEnd).thenComparing(EpochSamplingSkylineSegment::getSegmentType));


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Sanity checks - not necessary here, probably can be removed
        // First segment should start at 0
        if (segments.get(0).getStart() != 0)
            throw new IllegalArgumentException("First segment does not start at 0.");
        // Last segment should stop at tMRCA
        if (segments.get(segments.size()-1).getEnd() != tree.getRoot().getHeight())
            throw new IllegalArgumentException("Last segment does not end at tMRCA.");
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // Traverse times and types to get lineageCount, start and width
        int lineageCount = 0;
        double start = 0;
        for (EpochSamplingSkylineSegment segment : segments) {

            // Set segment
            segment.setStart(start);
            segment.setWidth(segment.getEnd()-start);
            segment.setLineageCount(lineageCount);

            // Update to next segment
            start = segment.getEnd();
            switch (segment.getSegmentType()) {
                case SAMPLE:
                    lineageCount++;
                    break;

                case COALESCENT:
                    lineageCount--;
                    break;

                default:
                    break;
            }
        }

        // Print the result
        // for (EpochSamplingSkylineSegment segment : segments) {
        //      System.out.println(segment);
        // }

        return segments;
    }


    // Use globals segments, epochtimes, groupsizes, popsizes, samplingintensities
    public boolean getSegments() {

        double currTime, prevTime = 0.0;


        // Collect the times and types of all segments
        segments = collectSegmentTimes();

        // Update arrays
        updateArrays();

        // Fill parameter values in each segment
        int groupIndex = 0,   // popSizes
            epochIndex = 0,   // samplingIntensities
            i = 0,  // Counts all segments
            j = 0;  // Counts only tree segments (
        for (EpochSamplingSkylineSegment segment : segments) {

            // popSize
            //if (j >= cumulativeGroupSizes[groupIndex])
            //    groupIndex++;
            if (i > groupLastSegments[groupIndex]) {
                currTime = segments.get(groupLastSegments[groupIndex]).getEnd();
                if (currTime - prevTime < minGroupWidth) {
                    Log.warning.println("Group "+groupIndex+" not wide enough (" + (currTime - prevTime) + ")");
                    return false;
                }
                prevTime = currTime;
                groupIndex++;
            }
            segment.setPopSize(popSizes.getValue(groupIndex));

            // samplingIntensity
            if (i > epochLastSegments[epochIndex])
                epochIndex++;
            //segment.setSamplingIntensity(i > epochFirstSamples[epochIndex] && i <= epochLastSamples[epochIndex] ? samplingIntensity.getValue(epochIndex) : 0.0);
            segment.setSamplingIntensity(i > epochLastSamples[epochIndex] ? 0.0 : samplingIntensity.getValue(epochIndex));

            if (segment.getSegmentType() == SegmentType.SAMPLE || segment.getSegmentType() == SegmentType.COALESCENT)
                j++;
            i++;

            //System.out.println(i+"\t"+j+"\t"+groupIndex+"\t"+epochIndex);
        }

        segmentsUpdated = true;

        return true;
    }


    @Override
    public double calculateLogP() {

        if (!segmentsUpdated) {
            if (!getSegments())
                return Double.NEGATIVE_INFINITY;
        }

        logP = 0.0;
        for (EpochSamplingSkylineSegment segment : segments) {
            logP += segment.calculateLikelihood();
        }
        return logP;
    }


    // Works better when outsourced to the SkylineSegment class
    /**
     * Calculates the log-likelihood of an interval under the coalescent with constant population size
     * (Copied and simplified from BayesianSkyline.java)
     *
     * @param popSize
     * @param beta
     * @param width
     * @param lineageCount
     * @param type
     * @return
     *
    public static double calculateIntervalLikelihood(double popSize, double beta, double width, int lineageCount, SegmentType type) {

        final double kchoose2 = Binomial.choose2(lineageCount);

        double lk = -width* (kchoose2/popSize + beta*popSize);

        switch (type) {
            case COALESCENT:
                lk += -Math.log(popSize);
                break;

            case SAMPLE:
                lk += Math.log(beta*popSize);
                break;

            default:
                break;
        }
        //System.out.println(lk + "-\t" + type+"\t"+lineageCount);

        return lk;
    }*/


    public String toString() {

        String outstr = String.format("%10s |%10s%10s%10s |%10s%12s |%10s%18s\n"+
                                      "-------------------------------------------------------------------------------------------------\n",
                                      "segment", "start", "end", "width", "lineages", "segmentType",
                                      "popSize", "samplingIntensity");

        if (!segmentsUpdated) {
            getSegments();
        }

        int i = 0;
        for (EpochSamplingSkylineSegment segment : segments) {
            outstr += String.format("%10d |",i) + segment.toString()+"\n";
            i++;
        }

        return outstr;
    }



    /************************************/
    /* Methods for getting change times */
    /************************************/

    public int getDimension() {
        return popSizes.getDimension();
    }

    /*
     * Return the i'th change time of the skyline parameter (for logging and skyline reconstruction)
     *
     * @param i
     * @return
     */
    public double getChangeTime(int i) {
        if (!segmentsUpdated) {
            getSegments();
        }
        return segments.get(groupLastSegments[i]).getEnd();
    }


    /**
     * TODO: Check boundary conditions
     *
     * @param t
     * @return
     *
    public double getPopSize(double t) {

        Collections.binarySearch(segments, t, (EpochSamplingSkylineSegment o1, double o2) -> Double.compare(o1.getEnd(),o2));

        segments.binarySearch(Comparator.comparing(EpochSamplingSkylineSegment::getEnd))
        int groupIndex = Arrays.binarySearch(   groupTimes, t);
        if (groupIndex < 0) {
            groupIndex = -groupIndex - 1;
        } else {
            groupIndex++;
        }

        return popSizes.getValue(Math.min(groupIndex, popSizes.getDimension()-1));
    }*/



    /****************************/
    /* Calculation Node methods */
    /****************************/

    @Override
    protected boolean requiresRecalculation() {
        segmentsUpdated = false;
        return true;
    }

    @Override
    public void store() {
        segmentsUpdated = false;
        super.store();
    }

    @Override
    public void restore() {
        segmentsUpdated = false;
        super.restore();
    }

}
