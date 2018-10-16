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
 * Bayesian Skyline plot with preferential sampling
 *
 * TODO: Check if prepared or not before calling prepare() to save time
 * TODO: Make it work with a tree input
 *
 * @author Louis du Plessis
 * @date 2018/09/18
 */
public class PreferentialBayesianSkyline extends TreeDistribution {

    final public Input<RealParameter> popSizeInput =
            new Input<>("popSizes","Effective population size (skyline parameter)");

    final public Input<RealParameter> samplingIntensityInput =
            new Input<>("samplingIntensity","Constant sampling intensity");

    final public Input<IntegerParameter> groupSizeInput =
            new Input<>("groupSizes", "The number of events in each group in the skyline (use robust design if not provided)");


    protected RealParameter popSizes, samplingIntensity;
    protected IntegerParameter groupSizes;
    protected int [] cumulativeGroupSizes;
    protected TreeIntervals intervals;
    protected double oldestSample;
    //protected double [] changeTimes;


    @Override
    public void initAndValidate() {

        int nrEvents, nrGroups;

        popSizes = popSizeInput.get();
        nrGroups = popSizes.getDimension();

        samplingIntensity = samplingIntensityInput.get();

        // Get tree intervals
        intervals = treeIntervalsInput.get();
        nrEvents  = intervals.getIntervalCount();

        if (groupSizeInput.get() == null) {
            // If groupSizes are not specified use robust design
            groupSizes = getRobustGroupSizes(nrEvents, nrGroups);
        } else {
            // Check groupSizes
            groupSizes = groupSizeInput.get();
        }

        // Number of groups does not equal the dimension of the skyline parameter
        if (groupSizes.getDimension() != nrGroups) {
            throw new IllegalArgumentException("Number of groups should match the dimension of the effective population size.");
        }

        // More groups than coalescent events
        if (groupSizes.getDimension() > nrEvents) {
            throw new IllegalArgumentException("There are more groups than events in the tree.");
        }

        cumulativeGroupSizes = new int[nrGroups];
        prepareTimes();

        // GroupSizes needs to add up to coalescent events
        if (cumulativeGroupSizes[nrGroups - 1] != nrEvents) {
            Log.warning.println("WARNING: The sum of the initial group sizes does not match the number of coalescent " +
                    "events in the tree. Initializing to equal group sizes");

            groupSizes.assignFromWithoutID(getRobustGroupSizes(nrEvents, nrGroups));

            // Recalculate cumulative group sizes, because group sizes have been changed
            prepareTimes();
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

        // Get oldest sample
        double currentTime = 0.0;
        oldestSample = 0.0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            currentTime += intervals.getInterval(i);

            if (intervals.getIntervalType(i) == IntervalType.SAMPLE) {
                oldestSample = currentTime;
            }
        }

    }

    @Override
    // TODO Replace eventIndex with i
    public double calculateLogP() {

        int    groupIndex = 0,
               eventIndex = 0;
        double width,
               currentTime = 0.0,
               currentPopSize,
               currentsamplingIntensity;
        // Get interval times

        prepareTimes();

        // Get likelihood for each segment
        logP = 0.0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {

            eventIndex++;
            if (eventIndex > cumulativeGroupSizes[groupIndex])
                groupIndex++;

            currentPopSize = popSizes.getArrayValue(groupIndex);
            currentsamplingIntensity = currentTime < oldestSample ? samplingIntensity.getValue() : 0.0;

            width = intervals.getInterval(i);

            //System.out.println(i+"\t"+groupIndex+"\t"+currentTime+"\t"+intervals.getInterval(i)+"\t"+currentPopSize+"\t"+currentsamplingIntensity);
            //System.out.println(i+"\t"+coalIndex+"\t"+groupIndex+"\t"+currentTime+"\t"+intervals.getInterval(i)+"\t"+intervals.getCoalescentEvents(i));

            // Calculate likelihood for interval
            logP += calculateIntervalLikelihood(currentPopSize, currentsamplingIntensity, width, intervals.getLineageCount(i), intervals.getIntervalType(i));

            currentTime += width;

        }

        return logP;
    }


    /**
     * Copied from BayesianSkyline.java (modified)
     *
     * @param popSize
     * @param width
     * @param lineageCount
     * @param type
     * @return
     */
    public static double calculateIntervalLikelihood(double popSize, double beta, double width, int lineageCount, IntervalType type) {

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

        return lk;
    }


    public int getDimension() {
        return popSizes.getDimension();
    }



    /* Calculation Node methods */

    @Override
    // TODO
    protected boolean requiresRecalculation() {
        return true;
    }

    public double getChangeTime(int i) {
        // check if prepared
        prepareTimes();

        int idx = Math.min(intervals.getIntervalCount()-1, cumulativeGroupSizes[i]);
        return intervals.getIntervalTime(idx);
    }


}
