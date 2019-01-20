package bsp.distributions;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.coalescent.IntervalType;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.math.Binomial;

import static beast.evolution.tree.coalescent.IntervalType.COALESCENT;



/**
 * The most basic Bayesian Skyline Plot implementation.
 *
 * No bells or whistles to facilitate anything. This implementation will return the same likelihood as
 * beast.evolution.tree.coalescent.BayesianSkyline, but will hopefully also be faster
 *
 * TODO: Write a version that doesn' use beast.evolution.tree.coalescent.TreeIntervals as this class is probably more
 * expensive than it needs to be.
 *
 * @author Louis du Plessis
 * @date 2019/01/16
 */
public class SimpleBSP extends TreeDistribution {


    final public Input<RealParameter> popSizeInput =
            new Input<>("popSizes","Effective population size (skyline parameter)");

    final public Input<IntegerParameter> groupSizeInput =
            new Input<>("groupSizes", "The number of events in each group in the skyline (use robust design if not provided)");

    protected TreeIntervals intervals;
    protected RealParameter popSizes;
    protected IntegerParameter groupSizes;

    protected int []    cumulativeGroupSizes,
            storedCumulativeGroupSizes;
    // protected double [] coalescentTimes, groupTimes;

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
            // Doesn't work... tree intervals are not updated with tree
            // intervals.initByName("tree", treeInput.get());
            throw new IllegalArgumentException("Only tree intervals (not tree) should be specified!");
        } else {
            intervals = treeIntervalsInput.get();
        }
        nrCoal = intervals.getSampleCount();


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
        //groupTimes[0]           = coalescentTimes[cumulativeGroupSizes[0]-1];
        for (int i = 1; i < cumulativeGroupSizes.length; i++) {
            cumulativeGroupSizes[i] = cumulativeGroupSizes[i-1] + groupSizes.getValue(i);
            //groupTimes[i]           = coalescentTimes[cumulativeGroupSizes[i]-1];
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
        for (int i = 0; i < intervals.getIntervalCount(); i++) {

            if (intervals.getIntervalType(i) == COALESCENT) {
                coalIndex++;
                if (coalIndex > cumulativeGroupSizes[groupIndex])
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
        arraysUpdated = false;
        return true;
    }

    @Override
    public void store() {
        arraysUpdated = false;
        //System.arraycopy(cumulativeGroupSizes, 0, storedCumulativeGroupSizes, 0, cumulativeGroupSizes.length);
        super.store();
    }

    @Override
    public void restore() {
        arraysUpdated = false;
        //int [] tmp = storedCumulativeGroupSizes;
        //storedCumulativeGroupSizes = cumulativeGroupSizes;
        //cumulativeGroupSizes = tmp;

        super.restore();
    }


}
