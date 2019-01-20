package bsp.distributions;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;
import beast.math.Binomial;
import beast.util.HeapSort;

/**
 * The most basic Bayesian Skyline Plot implementation.
 *
 * No bells or whistles to facilitate anything. This implementation will return the same likelihood as
 * beast.evolution.tree.coalescent.BayesianSkyline, but will hopefully also be faster
 *
 * Unlike SimpleBSP this version does not use the TreeIntervals class, but calculates intervals by itself
 *
 * @author Louis du Plessis
 * @date 2019/01/16
 */
public class SimpleBSP2  extends TreeDistribution {


    final public Input<RealParameter> popSizeInput =
            new Input<>("popSizes","Effective population size (skyline parameter)");

    final public Input<IntegerParameter> groupSizeInput =
            new Input<>("groupSizes", "The number of events in each group in the skyline (use robust design if not provided)");

    protected TreeInterface tree;
    protected RealParameter popSizes;
    protected IntegerParameter groupSizes;

    protected int []    cumulativeGroupSizes,
                  storedCumulativeGroupSizes;

    protected double      [] intervalTimes;
    protected double      [] intervalWidths;
    protected int         [] lineageCounts;
    protected SegmentType [] intervalTypes;

    protected double      [] storedIntervalTimes;
    protected double      [] storedIntervalWidths;
    protected int         [] storedLineageCounts;
    protected SegmentType [] storedIntervalTypes;

    protected boolean arraysUpdated = false;
    protected boolean intervalsKnown = false;

    @Override
    public void initAndValidate() {

        int nrEvents, nrCoal, nrGroups;

        //////////////////////////
        // Get skyline parameter
        popSizes = popSizeInput.get();
        nrGroups = popSizes.getDimension();


        ///////////////////////
        // Get tree intervals
        ///////////////////////
        // Get tree
        if (treeInput.get() == null) {
            throw new IllegalArgumentException("Only tree input (not tree intervals) should be specified!");
        } else {
            tree = treeInput.get();
        }
        nrCoal   = tree.getLeafNodeCount()-1;
        nrEvents = tree.getNodeCount();


        ////////////////////
        // Get group sizes
        // If groupSizes are not specified use robust design (equal group sizes)
        if (groupSizeInput.get() != null) {
            groupSizes = groupSizeInput.get();
        } else {
            groupSizes = getRobustGroupSizes(nrCoal, nrGroups);
        }

        // Group sizes does not equal the dimension of the skyline parameter
        if (groupSizes.getDimension() != nrGroups) {
            throw new IllegalArgumentException("Number of groups should match the dimension of the effective population size.");
        }

        // More groups than coalescent events
        if (groupSizes.getDimension() > nrCoal) {
            throw new IllegalArgumentException("There are more groups than coalescent events in the tree.");
        }

        cumulativeGroupSizes = new int[nrGroups];
        storedCumulativeGroupSizes = new int[nrGroups];
        // groupTimes = new double[nrGroups];

        intervalTimes  = new double[nrEvents];
        intervalWidths = new double[nrEvents];
        lineageCounts  = new int[nrEvents];
        intervalTypes  = new SegmentType[nrEvents];

        storedIntervalTimes  = new double[nrEvents];
        storedIntervalWidths = new double[nrEvents];
        storedLineageCounts  = new int[nrEvents];
        storedIntervalTypes  = new SegmentType[nrEvents];

        updateArrays();

        // GroupSizes needs to add up to coalescent events
        if (cumulativeGroupSizes[nrGroups - 1] != nrCoal) {
            Log.warning.println("WARNING: The sum of the initial group sizes does not match the number of coalescent " +
                    "events in the tree. Initializing to equal group sizes (robust design)");

            groupSizes.assignFromWithoutID(getRobustGroupSizes(nrCoal, nrGroups));

            // Recalculate cumulative group sizes, because group sizes have been changed
            updateArrays();
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
     * Updates the arrays used in likelihood calculation and other methods
     */
    protected void updateArrays() {


        // Get coalescent times (not necessary for likelihood computation)
        // coalescentTimes = intervals.getCoalescentTimes(coalescentTimes);

        // Get cumulative group sizes and times
        cumulativeGroupSizes[0] = groupSizes.getValue(0);
        for (int i = 1; i < cumulativeGroupSizes.length; i++) {
            cumulativeGroupSizes[i] = cumulativeGroupSizes[i-1] + groupSizes.getValue(i);
        }

        // Get the tree interval times, types and number of lineages per interval
        if (!intervalsKnown) {
            // 1) Get node order by time
            double[] times = new double[intervalTimes.length];
            SegmentType[] types = new SegmentType[intervalTimes.length];

            Node[] nodes = tree.getNodesAsArray();
            for (int i = 0; i < nodes.length; i++) {
                times[i] = nodes[i].getHeight();
                types[i] = nodes[i].isLeaf() ? SegmentType.SAMPLE : SegmentType.COALESCENT;
            }

            /**** SORTING ****/
            int[] indices = new int[intervalTimes.length];
            HeapSort.sort(times, indices);


            // 2) Traverse nodes from present to TMRCA while reordering and also calculate lineageCounts and intervalWidths
            int lineageCount = 1;
            intervalTimes[0] = times[indices[0]];
            intervalTypes[0] = types[indices[0]];
            intervalWidths[0] = 0;
            lineageCounts[0] = 0;
            for (int i = 1; i < intervalTimes.length; i++) {

                intervalTimes[i] = times[indices[i]];
                intervalTypes[i] = types[indices[i]];
                intervalWidths[i] = intervalTimes[i] - intervalTimes[i - 1];
                lineageCounts[i] = lineageCount;

                switch (intervalTypes[i]) {
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
            intervalsKnown = true;
        }

        arraysUpdated = true;
    }

    @Override
    public double calculateLogP() {

        int    groupIndex = 0,
                coalIndex  = 0;
        double currentPopSize;

        // Update arrays
        if (!arraysUpdated) {
            updateArrays();
        }

        // Get likelihood for each segment
        logP = 0.0;
        for (int i = 0; i < intervalTypes.length; i++) {

            if (intervalTypes[i] == SegmentType.COALESCENT) {
                coalIndex++;
                if (coalIndex > cumulativeGroupSizes[groupIndex])
                    groupIndex++;
            }

            currentPopSize = popSizes.getArrayValue(groupIndex);
            logP += calculateIntervalLikelihood(currentPopSize, intervalWidths[i], lineageCounts[i], intervalTypes[i]);
        }

        return logP;
    }


    /**
     * Calculates the log-likelihood of an interval under the coalescent with constant population size
     * (Copied and simplified from BayesianSkyline.java)
     *
     * @param popSize
     * @param width
     * @param lineageCount
     * @param type
     * @return
     */
    public static double calculateIntervalLikelihood(double popSize, double width, int lineageCount, SegmentType type) {

        final double kchoose2 = Binomial.choose2(lineageCount);

        double lk = -kchoose2 * width / popSize;

        switch (type) {
            case COALESCENT:
                lk += -Math.log(popSize);

                break;
            default:
                break;
        }

        //System.out.printf("%10.5f %10s |%10d%12s |%10s\n",
        //        lk, width, lineageCount, type, popSize);


        return lk;
    }


    /************************************/
    /* Methods for getting change times */
    /************************************/

    public int getDimension() {
        return popSizes.getDimension();
    }

    /**
     * Return the i'th change time of the skyline parameter (for logging and skyline reconstruction)
     *
     * @param i
     * @return
     *
    public double getChangeTime(int i) {
    if (!arraysUpdated) {
    updateArrays();
    }

    return groupTimes[i];
    }*/


    /**
     * TODO: Check boundary conditions
     *
     * @param t
     * @return
     *
    public double getPopSize(double t) {

    int groupIndex = Arrays.binarySearch(groupTimes, t);
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
        arraysUpdated  = false;
        intervalsKnown = !tree.somethingIsDirty();

        return true;
    }

    @Override
    public void store() {
        //arraysUpdated = false;

        System.arraycopy(cumulativeGroupSizes, 0, storedCumulativeGroupSizes, 0, cumulativeGroupSizes.length);
        System.arraycopy(intervalTimes, 0, storedIntervalTimes, 0, intervalTimes.length);
        System.arraycopy(intervalTypes, 0, storedIntervalTypes, 0, intervalTypes.length);
        System.arraycopy(intervalWidths, 0, storedIntervalWidths, 0, intervalWidths.length);
        System.arraycopy(lineageCounts, 0, storedLineageCounts, 0, lineageCounts.length);

        super.store();
    }

    @Override
    public void restore() {
        //arraysUpdated = false;

        int [] tmp = storedCumulativeGroupSizes;
        storedCumulativeGroupSizes = cumulativeGroupSizes;
        cumulativeGroupSizes = tmp;

        double [] tmp2 = storedIntervalTimes;
        storedIntervalTimes = intervalTimes;
        intervalTimes = tmp2;

        double [] tmp3 = storedIntervalWidths;
        storedIntervalWidths = intervalWidths;
        intervalWidths = tmp3;

        int [] tmp4 = storedLineageCounts;
        storedLineageCounts = lineageCounts;
        lineageCounts = tmp4;

        SegmentType [] tmp5 = storedIntervalTypes;
        storedIntervalTypes = intervalTypes;
        intervalTypes = tmp5;

        super.restore();
    }

}
