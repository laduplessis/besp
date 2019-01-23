package bsp.distributions;


import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.coalescent.IntervalType;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.math.Binomial;
import beast.util.HeapSort;
import beast.util.Randomizer;

import java.util.Arrays;

import static beast.evolution.tree.coalescent.IntervalType.COALESCENT;

/**
 * Bayesian Skyline Plot implementation.
 *
 * This implementation will return the same likelihood as
 * beast.evolution.tree.coalescent.BayesianSkyline, but will hopefully also be faster
 *
 * This implementation makes use of the TreeIntervals class, like the original BayesianSkyline
 *
 * Only one parameter (popSizes size) which can only change at coalescent times, which are grouped with
 * popSizeGroupSizes.
 * Times for groups are logged as well and there is a minimum width for each group.
 *
 * @author Louis du Plessis
 * @date 2019/01/20
 */
public class BSP extends TreeDistribution {

    final public Input<RealParameter> popSizeInput =
            new Input<>("popSizes","Effective population size (skyline parameter)");

    final public Input<IntegerParameter> popSizeGroupSizeInput =
            new Input<>("popSizeGroupSizes", "The number of events in each population size group in the skyline (use robust design if not provided)");

    final public Input<Double> minWidthInput =
            new Input<>("minWidth","Minimum width of a group (end-start)",0.0);

    final public Input<Integer> numInitializationAttempsInput =
            new Input<>("numInitializationAttempts","Number of times to try to initialize the group sizes"+
                              "if the minimum group width constraint is not satisfied",10000);

    protected TreeIntervals intervals;
    protected RealParameter popSizes;
    protected IntegerParameter popSizeGroupSizes;

    protected int []    cumulativePopSizeGroupSizes;
                //storedCumulativePopSizeGroupSizes;
    protected double [] popSizeGroupTimes;
    protected double minWidth;

    protected boolean arraysUpdated = false;


    @Override
    public void initAndValidate() {

        int nrCoal, nrGroups;

        //////////////////////////
        // Get skyline parameter
        popSizes = popSizeInput.get();
        nrGroups = popSizes.getDimension();


        ///////////////////////
        // Get tree intervals
        if (treeInput.get() != null) {
            throw new IllegalArgumentException("Only tree intervals (not tree) should be specified!");
        } else {
            intervals = treeIntervalsInput.get();
        }
        nrCoal = intervals.getSampleCount();


        ////////////////////////////
        // Get minimum group time
        minWidth = minWidthInput.get();


        ////////////////////
        // Get group sizes
        // If popSizeGroupSizes are not specified use robust design (equal group sizes)
        if (popSizeGroupSizeInput.get() != null) {
            popSizeGroupSizes = popSizeGroupSizeInput.get();
        } else {
            popSizeGroupSizes = getRobustpopSizeGroupSizes(nrCoal, nrGroups, 1, Integer.MAX_VALUE);
        }

        // Group sizes does not equal the dimension of the skyline parameter
        if (popSizeGroupSizes.getDimension() != nrGroups) {
            throw new IllegalArgumentException("Number of groups should match the dimension of the effective population size.");
        }

        // More groups than coalescent events
        if (popSizeGroupSizes.getDimension() > nrCoal) {
            throw new IllegalArgumentException("There are more groups than coalescent events in the tree.");
        }

        /////////////////////
        // Initialise arrays
        cumulativePopSizeGroupSizes = new int[nrGroups];
        //storedCumulativepopSizeGroupSizes = new int[nrGroups];
        popSizeGroupTimes = new double[nrGroups];
        updateArrays();

        // popSizeGroupSizes needs to add up to coalescent events
        if (cumulativePopSizeGroupSizes[nrGroups - 1] != nrCoal) {
            Log.warning.println("WARNING: The sum of the initial group sizes does not match the number of coalescent " +
                    "events in the tree. Initializing to equal group sizes (robust design)");

            popSizeGroupSizes.assignFromWithoutID(getRobustpopSizeGroupSizes(nrCoal, nrGroups,
                                                                             popSizeGroupSizes.getLower(), popSizeGroupSizes.getUpper()));

            // Recalculate cumulative group sizes, because group sizes have been changed
            updateArrays();
        }

        int i = 0;
        int numInitializationAttemps = numInitializationAttempsInput.get();
        while (!checkGroupWidths(popSizeGroupTimes, minWidth)) {
            if (i > numInitializationAttemps) {
                throw new IllegalArgumentException("Minimum group width is still shorter than minWidth (" + minWidth + ") "
                        + "after "+numInitializationAttemps+" attempts at redistributing group sizes.\n"
                        + "Try decreasing the number of groups or the minimum group width.");
                        //+ "Current group times: " + Arrays.toString(popSizeGroupTimes));
            }

            redistributeGroups(popSizeGroupSizes, popSizeGroupTimes);
            updateArrays();
            i++;
        }

    }

    /**
     * Updates the arrays used in likelihood calculation and other methods
     */
    protected void updateArrays() {

        // Get coalescent times (not necessary for likelihood computation)
        double [] coalescentTimes = intervals.getCoalescentTimes(null);

        // Get cumulative group sizes and times
        cumulativePopSizeGroupSizes[0] = popSizeGroupSizes.getValue(0);
        popSizeGroupTimes[0]           = coalescentTimes[cumulativePopSizeGroupSizes[0]-1];
        for (int i = 1; i < cumulativePopSizeGroupSizes.length; i++) {
            cumulativePopSizeGroupSizes[i] = cumulativePopSizeGroupSizes[i-1] + popSizeGroupSizes.getValue(i);
            popSizeGroupTimes[i]           = coalescentTimes[cumulativePopSizeGroupSizes[i]-1];
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

            if (!checkGroupWidths(popSizeGroupTimes, minWidth)) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        // Get likelihood for each segment
        logP = 0.0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {

            if (intervals.getIntervalType(i) == COALESCENT) {
                coalIndex++;
                if (coalIndex > cumulativePopSizeGroupSizes[groupIndex])
                    groupIndex++;
            }

            currentPopSize = popSizes.getArrayValue(groupIndex);
            logP += calculateIntervalLikelihood(currentPopSize, intervals.getInterval(i), intervals.getLineageCount(i), intervals.getIntervalType(i));
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
    public static double calculateIntervalLikelihood(double popSize, double width, int lineageCount, IntervalType type) {

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

    @Override
    public String toString() {

        String outstr = this.getID();



        return outstr;
    }


    /*******************************/
    /* Methods for population size */
    /*******************************/

    public int getPopSizeDimension() {
        return popSizes.getDimension();
    }

    /**
     * Return the i'th change time of the skyline parameter (for logging and skyline reconstruction)
     *
     * @param i
     * @return
     */
    public double getPopSizeChangeTime(int i) {
        if (!arraysUpdated) {
            updateArrays();
        }

        return popSizeGroupTimes[i];
    }


    /**
     * TODO: Check boundary conditions
     *
     * @param t
     * @return
     */
    public double getPopSize(double t) {

        int groupIndex = Arrays.binarySearch(popSizeGroupTimes, t);
        if (groupIndex < 0) {
            groupIndex = -groupIndex - 1;
        } else {
            groupIndex++;
        }

        return popSizes.getValue(Math.min(groupIndex, popSizes.getDimension()-1));
    }


    /****************************/
    /* Calculation Node methods */
    /****************************/

    @Override
    protected boolean requiresRecalculation() {
        arraysUpdated = false;
        return true;
    }

    @Override
    public void store() {
        arraysUpdated = false;
        //System.arraycopy(cumulativePopSizeGroupSizes, 0, storedCumulativepopSizeGroupSizes, 0, cumulativePopSizeGroupSizes.length);
        super.store();
    }

    @Override
    /**
     * No real speed advantage to storing and restoring arrays instead of just always updating
     */
    public void restore() {
        arraysUpdated = false;
        //int [] tmp = storedCumulativepopSizeGroupSizes;
        //storedCumulativepopSizeGroupSizes = cumulativePopSizeGroupSizes;
        //cumulativePopSizeGroupSizes = tmp;

        super.restore();
    }



    /*****************************************************/
    /* Methods for checking and initialising group sizes */
    /*****************************************************/




    protected boolean checkGroupWidths(double [] groupTimes, double minWidth) {

        double width, prev = 0.0;

        if (minWidth > 0.0) {
            for (int i = 0; i < groupTimes.length; i++) {
                width = groupTimes[i]-prev;
                prev  = groupTimes[i];
                if (width < minWidth) {
                    // System.out.println("Duration of group " + (i+1) + " too short: " + width);
                    return false;
                }

            }
        }
        return true;
    }


    /**
     * This method should not be necessary as no operator should adjust groups outside of their limits
     *
     * @return
     */
    protected boolean checkGroupSizes(IntegerParameter groupSizes) {

        if (groupSizes.getLower() > 0 || groupSizes.getUpper() < Integer.MAX_VALUE-1) {
            for (Integer group : groupSizes.getValues()) {
                if (group < groupSizes.getLower() || group > groupSizes.getUpper()) {
                    return false;
                }
            }
        }
        return true;

    }



    /**
     * Distribute an equal number of events in each group
     *
     * @param events
     * @param groups
     * @return
     */
    protected IntegerParameter getRobustpopSizeGroupSizes(int events, int groups, int lower, int upper) {

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
        parameter.setBounds(Math.max(1,lower), Math.min(Integer.MAX_VALUE,upper));

        return(parameter);
    }


    /**
     * (works similar to a DeltaExchangeOperator
     *
     * Select two random indices of sizes (two randomly selected groups)
     * Then select random delta
     * Increase one group size by delta, decrease the other by delta
     *
     * Delta is adjusted so the parameter respects the limits
     * All groups > 0
     * No group >= sum(groups)
     *
     * @param sizes
     */
    protected Integer [] redistributeGroupsRandom(Integer [] sizes) {

        int minSize = 0, maxSize = 0;
        for (int i = 0; i < sizes.length; i++) {
            maxSize += sizes[i];
        }

        int dim1 = Randomizer.nextInt(sizes.length);
        int dim2 = Randomizer.nextInt(sizes.length-1);
        if (dim2 >= dim1) {
            ++dim2;
        }

        // Prevents adjusted sizes to be < minSize or > maxSize
        int deltaMax = Math.min(maxSize - sizes[dim1], sizes[dim2] - minSize) - 1;

        if (deltaMax > 0) {
            int delta = Randomizer.nextInt((int) Math.round(deltaMax)) + 1;

            sizes[dim1] += delta;
            sizes[dim2] -= delta;

            if (sizes[dim1] <= minSize || sizes[dim1] >= maxSize || sizes[dim2] <= minSize || sizes[dim2] >= maxSize) {
                throw new RuntimeException("Problem with distributing groups. This should never happen.");
            }
        }
        return sizes;
    }

    /**
     * Draw an integer from [0, weights.length) with associated weights.
     *
     * @param weights
     * @return
     */
    private int weightedRandomSample(double [] weights) {

        // Turn weights into cumulative weights
        for (int i = 1; i < weights.length; i++) {
            weights[i] += weights[i-1];
        }

        // Normalise so the sum is equal to 1
        for (int i = 0; i < weights.length; i++) {
            weights[i] /= weights[weights.length-1];
        }

        // Sample random nr and find where it fits
        double r = Randomizer.nextDouble();
        return Math.abs(Arrays.binarySearch(weights, r))-1;
    }


    /**
     * Heuristic random algorithm for redistributing groups to maximise the minimum group duration
     *
     * 1) Sort groups by durations
     * 2) To select dim1, always select the group with the shortest duration
     * 3) To select dim2, draw weighted random sample, weights are duration where all groups of size 1 have weight 0
     * 4) Always set delta to 1 (alternatively, select delta randomly between 1 and max possible delta)
     * 5) sizes[dim1] += delta, sizes[dim2] -= delta
     *
     */
    protected void redistributeGroups(IntegerParameter groupSizes, double [] groupTimes) {

        int dim1, dim2;

        // Get group widths
        double[] widths = new double[groupTimes.length];
        double prev = 0.0;
        for (int i = 0; i < widths.length; i++) {
            widths[i] = groupTimes[i] - prev;
            prev      = groupTimes[i];
        }

        // Sort group widths and get order
        int[] indices = new int[widths.length];
        HeapSort.sort(widths, indices);

        // Set dim1 to group with the shortest width, draw dim2 weighted by the width
        Integer[] sizes = groupSizes.getValues();
        dim1 = indices[0];

        double [] weights = new double [widths.length];
        for (int i = 0; i < widths.length; i++) {
            weights[i] = sizes[indices[i]] > groupSizes.getLower() ? widths[indices[i]] : 0;
        }
        dim2 = indices[weightedRandomSample(weights)];


        // Prevents adjusted sizes to be < minSize or > maxSize
        int delta,
            deltaMax,
            minSize = 0,
            maxSize = 0;
        for (int i = 0; i < sizes.length; i++) {
            maxSize += sizes[i];
        }
        deltaMax = Math.min(maxSize - sizes[dim1], sizes[dim2] - minSize) - 1;

        if (deltaMax > 0) {
            //delta = Randomizer.nextInt((int) Math.round(deltaMax)) + 1;
            delta = 1;

            sizes[dim1] += delta;
            sizes[dim2] -= delta;

            if (sizes[dim1] <= minSize || sizes[dim1] >= maxSize || sizes[dim2] <= minSize || sizes[dim2] >= maxSize) {
                throw new RuntimeException("Problem with limits while redistributing groups. This should never happen.");
            }

            for (int j = 0; j < sizes.length; j++) {
                groupSizes.setValue(j, sizes[j]);
            }
        }
    }


    /*
    protected void redistributeGroups() {


        double minTime = Double.POSITIVE_INFINITY,
                maxTime = 0,
                time;
        int dim1 = 0, dim2;

        Integer [] sizes = popSizeGroupSizes.getValues();



        dim2 = Randomizer.nextInt(popSizeGroupTimes.length-1);
        if (dim2 >= dim1) {
            ++dim2;
        }

        for (int i = 1; i < popSizeGroupTimes.length; i++) {
            time = popSizeGroupTimes[i] - popSizeGroupTimes[i-1];

            if (time < minTime) {
                minTime = time;
                dim1 = i;
            }

            if (time > maxTime && sizes[i] > 1) {
                maxTime = time;
                dim2 = i;
            }
        }

        System.out.println(dim1+"\t("+minTime+")\t\t"+dim2+"\t("+maxTime+")");

        // Prevents adjusted sizes to be < minSize or > maxSize
        int delta,
                deltaMax,
                minSize = 0,
                maxSize = 0;
        for (int i = 0; i < sizes.length; i++) {
            maxSize += sizes[i];
        }
        deltaMax = Math.min(maxSize - sizes[dim1], sizes[dim2] - minSize) - 1;

        //System.out.println(sizes[dim1] + "\t" + sizes[dim2] + "\t" + deltaMax);

        if (deltaMax > 0) {
            delta = Randomizer.nextInt((int) Math.round(deltaMax)) + 1;
            //delta = 1;

            sizes[dim1] += delta;
            sizes[dim2] -= delta;
            // System.out.println(Arrays.toString(sizes) + "\t" + delta + "\t" + deltaMax);

            if (sizes[dim1] <= minSize || sizes[dim1] >= maxSize || sizes[dim2] <= minSize || sizes[dim2] >= maxSize) {
                throw new RuntimeException("This should not happen!");
            }

            for (int j = 0; j < sizes.length; j++) {
                popSizeGroupSizes.setValue(j, sizes[j]);
            }

        } else {
            System.out.println("Cannot alter");
        }

    }
    */



}