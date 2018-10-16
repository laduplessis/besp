package bsp.distributions;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.coalescent.IntervalType;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.math.Binomial;

import java.util.Arrays;

/**
 *
 *
 * @author Louis du Plessis
 * @date 2018/09/11
 */
public class GeneralBayesianSkyline extends TreeDistribution {

    final public Input<RealParameter> popSizeInput =
            new Input<>("popSizes","Effective population size (skyline parameter)");

    final public Input<RealParameter> logPopSizeInput =
            new Input<>("logPopSizes","Log effective population size (skyline parameter)", Input.Validate.XOR, popSizeInput);

    final public Input<IntegerParameter> groupSizeInput =
            new Input<>("groupSizes", "The number of events in each group in the skyline (use robust design if not provided)");


    protected boolean useLog;         // Calculate logNe instead of Ne
    protected RealParameter popSizes, logPopSizes;
    protected IntegerParameter groupSizes;
    protected int [] cumulativeGroupSizes;
    protected double [] coalescentTimes; // actually unnecessary

    protected TreeIntervals intervals;
    // intervals.intervals[i] = width of interval[i]
    // intervals.lineageCounts[i] = nr lineages in interval[i]
    // intervals.coalescentTimes[i] = times of coalescent events
    //
    // intervals.getIntervalCount() intervals.getSampleCount()
    // intervals.getInterval(i)
    // intervals.getIntervalTime(i) = start of interval i
    // intervals.getIntervals() = return intervals
    // intervals.getIntervalType(i) = [ COALESCENT | SAMPLE | NOTHING ]
    // intervals.getCoalescentTimes() = times of coalescent events


    @Override
    public void initAndValidate() {

        int nrCoal, nrGroups;

        // Get skyline parameter
        if (popSizeInput.get() == null) {
            useLog = true;
            logPopSizes = logPopSizeInput.get();
            nrGroups    = logPopSizes.getDimension();
        } else {
            useLog = false;
            popSizes = popSizeInput.get();
            nrGroups = popSizes.getDimension();
        }

        // Get tree intervals
        intervals = treeIntervalsInput.get();
        nrCoal    = intervals.getSampleCount();

        if (groupSizeInput.get() == null) {

            // If groupSizes are not specified use robust design (not yet implemented)
            groupSizes = getRobustGroupSizes(nrCoal, nrGroups);

        } else {
            // Check groupSizes
            groupSizes = groupSizeInput.get();
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
        prepareTimes();

        // GroupSizes needs to add up to coalescent events
        if (cumulativeGroupSizes[nrGroups - 1] != nrCoal) {
            Log.warning.println("WARNING: The sum of the initial group sizes does not match the number of coalescent " +
                                "events in the tree. Initializing to equal group sizes");

            groupSizes.assignFromWithoutID(getRobustGroupSizes(nrCoal, nrGroups));
            //groupSizes = getRobustGroupSizes(nrCoal, nrGroups);
        }


    }

    /**
     *
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


    protected void prepareTimes() {

        // Get cumulative group sizes
        cumulativeGroupSizes[0] = groupSizes.getValue(0);
        for (int i = 1; i < cumulativeGroupSizes.length; i++) {
            cumulativeGroupSizes[i] = cumulativeGroupSizes[i-1] + groupSizes.getValue(i);
        }

        // Get coalescent times
        coalescentTimes = intervals.getCoalescentTimes(coalescentTimes);

        //



    }

    @Override
    public double calculateLogP() {

        int    groupIndex = 0,
               coalIndex = 0;
        double currentPopSize, currentLogPopSize,
               currentTime = 0.0;

        // Get interval times

        prepareTimes();

        // Get likelihood for each segment
        logP = 0.0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {

            // Set current population size to population size in the middle of the current interval
            //double currentPopSize2 = getPopSize(currentTime + intervals.getInterval(i)/2);

            if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                coalIndex++;
                if (coalIndex > cumulativeGroupSizes[groupIndex])
                    groupIndex++;
            }

            if (useLog) {
                //System.out.println(groupIndex);
                currentLogPopSize = logPopSizes.getArrayValue(groupIndex);
                logP += calculateIntervalLikelihoodLogPopSize(currentLogPopSize, intervals.getInterval(i), intervals.getLineageCount(i), intervals.getIntervalType(i));

            } else {
                currentPopSize = popSizes.getArrayValue(groupIndex);

                //System.out.println(i+"\t"+coalIndex+"\t"+groupIndex+"\t"+currentTime+"\t"+intervals.getInterval(i)+"\t"+currentPopSize+"\t"+currentPopSize2);

                //System.out.println(i+"\t"+coalIndex+"\t"+groupIndex+"\t"+currentTime+"\t"+intervals.getInterval(i)+"\t"+intervals.getCoalescentEvents(i));

                // Calculate likelihood for interval
                logP += calculateIntervalLikelihood(currentPopSize, intervals.getInterval(i), intervals.getLineageCount(i), intervals.getIntervalType(i));
            }

            currentTime += intervals.getInterval(i);
        }


        return logP;
    }


    /**
     * @param t time
     * @return
     *
    public double getPopSize(double t) {


        if (t > coalescentTimes[coalescentTimes.length - 1])
            return popSizes.getArrayValue(popSizes.getDimension() - 1);

        int epoch = Arrays.binarySearch(coalescentTimes, t);
        if (epoch < 0) {
            epoch = -epoch - 1;
        }

        int groupIndex = Arrays.binarySearch(cumulativeGroupSizes, epoch);

        if (groupIndex < 0) {
            groupIndex = -groupIndex - 1;
        } else {
            groupIndex++;
        }
        if (groupIndex >= popSizes.getDimension()) {
            groupIndex = popSizes.getDimension() - 1;
        }

        return popSizes.getArrayValue(groupIndex);
    }*/


    /**
     * Copied from BayesianSkyline.java (modified)
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

        return lk;
    }

    /**
     * Copied from BayesianSkyline.java (modified)
     *
     * @param logPopSize
     * @param width
     * @param lineageCount
     * @param type
     * @return
     */
    public static double calculateIntervalLikelihoodLogPopSize(double logPopSize, double width, int lineageCount, IntervalType type) {

        final double kchoose2 = Binomial.choose2(lineageCount);

        double lk = -kchoose2 * width * Math.exp(-logPopSize);

        switch (type) {
            case COALESCENT:
                lk += -logPopSize;

                break;
            default:
                break;
        }

        return lk;
    }


    /* Calculation Node methods */

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

}
