package bsp.distributions;


import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * A Bayesian Skyline Plot implementation using a list of segments for the tree intervals
 *
 *
 * @author Louis du Plessis
 * @date 2019/01/17
 */
public class SegmentedBSP extends TreeDistribution {


    final public Input<RealParameter> popSizeInput =
            new Input<>("popSizes","Effective population size (skyline parameter)", Input.Validate.REQUIRED);

    final public Input<IntegerParameter> groupSizeInput =
            new Input<>("groupSizes", "The number of events in each group in the skyline (use robust design if not provided)");


    protected TreeInterface tree;
    protected RealParameter popSizes;
    protected IntegerParameter groupSizes;

    protected int [] cumulativeGroupSizes,
                     groupLastSegments;

    protected int [] storedCumulativeGroupSizes,
                     storedGroupLastSegments;

    boolean segmentsUpdated = false,
            treeChanged     = true;

    // Segments of the skyline model
    List<BSPSegment>  segments,
                storedSegments;


    @Override
    public void initAndValidate() {

        int nrEvents, nrGroups;

        //////////////////////////
        // Get skyline parameter
        popSizes = popSizeInput.get();
        nrGroups = popSizes.getDimension();

        ///////////////////////
        // Get tree
        if (treeInput.get() == null) {
            throw new IllegalArgumentException("Only tree input (not tree intervals) should be specified!");
        } else {
            tree = treeInput.get();
        }
        nrEvents = tree.getLeafNodeCount()-1;

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

        storedCumulativeGroupSizes = new int[nrGroups];
        storedGroupLastSegments = new int [nrGroups];

        // GroupSizes needs to add up to coalescent events
        cumulativeGroupSizes[0] = groupSizes.getValue(0);
        for (int i = 1; i < cumulativeGroupSizes.length; i++) {
            cumulativeGroupSizes[i] = cumulativeGroupSizes[i-1] + groupSizes.getValue(i);
        }

        if (cumulativeGroupSizes[nrGroups - 1] != nrEvents) {
            Log.warning.println("WARNING: The sum of the initial group sizes does not match the number of coalescent " +
                    "events in the tree. Initializing to equal group sizes (robust design)");

            groupSizes.assignFromWithoutID(getRobustGroupSizes(nrEvents, nrGroups));

            // Recalculate cumulative group sizes, because group sizes have been changed
            //updateArrays();
        }


        treeChanged = true;
        getSegments();
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
     * Updates the arrays used in likelihood calculation and other methods
     *
     * - cumulateGroupSizes: cumulative sum of group sizes
     *
     */
    protected void updateArrays() {

        // Get cumulative group sizes
        cumulativeGroupSizes[0] = groupSizes.getValue(0);
        for (int i = 1; i < cumulativeGroupSizes.length; i++) {
            cumulativeGroupSizes[i] = cumulativeGroupSizes[i-1] + groupSizes.getValue(i);
        }
    }


    /**
     * 1.) Collect the end times and types of all segments
     * 2.) Orders times and types
     * 3.) Count lineages and widths of each segment
     *
     */
    protected List collectSegmentTimes() {

        Node[] nodes = tree.getNodesAsArray();
        ArrayList<BSPSegment> segments = new ArrayList<>();

        // Get coalescent and sampling events on the tree, epoch times and sort by end time, then segment type
        // such that for multiple segments ending at the same time, the epoch change time is the last segment
        for (int i = 0; i < nodes.length; i++) {
            segments.add(new BSPSegment(nodes[i].getHeight(),
                    nodes[i].isLeaf() ? SegmentType.SAMPLE : SegmentType.COALESCENT));
        }
        //Collections.sort(segments, (o1, o2) -> (Double.compare(o1.getEnd(),o2.getEnd())));
        segments.sort(Comparator.comparing(BSPSegment::getEnd).thenComparing(BSPSegment::getSegmentType));


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
        for (BSPSegment segment : segments) {

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


    // Use globals segments, groupsizes, popsizes
    public void getSegments() {

        //System.out.println("getSegments() "+treeChanged);

        // Collect the times and types of all segments
        if (treeChanged) segments = collectSegmentTimes();

        // Update arrays
        updateArrays();

        // Fill parameter values in each segment
        int groupIndex = 0,   // popSizes
            i = 0,            // Counts all segments
            j = 0;            // Counts coalescent events
        for (BSPSegment segment : segments) {

            segment.setPopSize(popSizes.getValue(groupIndex));

            if (segment.getSegmentType() == SegmentType.COALESCENT) {
                j++;
                if (j >= cumulativeGroupSizes[groupIndex]) {
                    groupLastSegments[groupIndex] = i;
                    groupIndex++;
                }
            }
            i++;
        }

        segmentsUpdated = true;
    }



    @Override
    public double calculateLogP() {

        if (!segmentsUpdated) {
            getSegments();
        }

        logP = 0.0;
        for (BSPSegment segment : segments) {
            logP += segment.calculateLikelihood();
        }
        return logP;
    }



    public String toString() {

        String outstr = String.format("%10s |%10s%10s%10s |%10s%12s\n"+
                        "------------------------------------------------------------------------------\n",
                "segment", "start", "end", "width", "lineages", "segmentType",
                "popSize");

        if (!segmentsUpdated) {
            getSegments();
        }

        int i = 0;
        for (BSPSegment segment : segments) {
            outstr += String.format("%10d |",i) + segment.toString()+"\n";
            i++;
        }

        outstr += "\nGroup sizes           : "+ Arrays.toString(groupSizes.getValues());
        outstr += "\nCumulative group sizes: "+ Arrays.toString(cumulativeGroupSizes);
        outstr += "\nGroup last segments   : "+ Arrays.toString(groupLastSegments)+"\n";

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



    /****************************/
    /* Calculation Node methods */
    /****************************/

    @Override
    protected boolean requiresRecalculation() {
        // Always update cumulative group sizes and the population sizes of segments
        // Only update the segment times, types and lineages etc. if the tree has changed
        segmentsUpdated = false;
        treeChanged = tree.somethingIsDirty();

        return true;
    }

    @Override
    public void store() {
        //segmentsUpdated = false;
        //treeChanged = true;

        System.arraycopy(cumulativeGroupSizes, 0, storedCumulativeGroupSizes, 0, cumulativeGroupSizes.length);
        System.arraycopy(groupLastSegments,    0, storedGroupLastSegments,    0, groupLastSegments.length);

        storedSegments = new ArrayList();
        for (BSPSegment segment : segments) {
            storedSegments.add(new BSPSegment(segment));
        }

        super.store();
    }

    @Override
    public void restore() {
        //segmentsUpdated = false;
        //treeChanged = true;

        int [] tmp = storedCumulativeGroupSizes;
        storedCumulativeGroupSizes = cumulativeGroupSizes;
        cumulativeGroupSizes = tmp;

        tmp = storedGroupLastSegments;
        storedGroupLastSegments = groupLastSegments;
        groupLastSegments = tmp;

        List<BSPSegment> tmp2 = storedSegments;
        storedSegments = segments;
        segments = tmp2;

        super.restore();
    }
}
