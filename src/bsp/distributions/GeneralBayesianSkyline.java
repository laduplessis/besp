package bsp.distributions;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.coalescent.IntervalType;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.math.Binomial;

import java.util.Arrays;

import static beast.evolution.tree.coalescent.IntervalType.COALESCENT;

/**
 *  Base Bayesian Skyline implementation
 *
 * @author Louis du Plessis
 * @date 2018/09/11
 */
@Description("Base implementation of the Bayesian skyline plot.")
@Citation(value="Drummond, A. J., Rambaut, A., Shapiro, B., & Pybus, O. G. (2005).\n" +
                "Bayesian coalescent inference of past population dynamics from molecular sequences.\n" +
                "Molecular biology and evolution, 22(5), 1185-1192.",
                year = 2005, firstAuthorSurname = "Drummond", DOI="10.1093/molbev/msi103")
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
    protected double [] coalescentTimes;

    protected TreeIntervals intervals;



    @Override
    public void initAndValidate() {

        int nrCoal, nrGroups;

        //////////////////////////
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

        ///////////////////////
        // Get tree intervals
        // Allow input tree or tree intervals to be specified. If the input TreeInterface is not a Tree it will crash.
        if (treeInput.get() != null) {
            intervals = new TreeIntervals((Tree) treeInput.get());
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
        prepareTimes();

        // GroupSizes needs to add up to coalescent events
        if (cumulativeGroupSizes[nrGroups - 1] != nrCoal) {
            Log.warning.println("WARNING: The sum of the initial group sizes does not match the number of coalescent " +
                                "events in the tree. Initializing to equal group sizes (robust design)");

            groupSizes.assignFromWithoutID(getRobustGroupSizes(nrCoal, nrGroups));

            // Recalculate cumulative group sizes, because group sizes have been changed
            prepareTimes();
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


    protected void prepareTimes() {

        // Get cumulative group sizes
        cumulativeGroupSizes[0] = groupSizes.getValue(0);
        for (int i = 1; i < cumulativeGroupSizes.length; i++) {
            cumulativeGroupSizes[i] = cumulativeGroupSizes[i-1] + groupSizes.getValue(i);
        }

        // Get coalescent times (unnecessary)
        coalescentTimes = intervals.getCoalescentTimes(coalescentTimes);

    }

    @Override
    public double calculateLogP() {

        int    groupIndex = 0,
               coalIndex  = 0;
        double currentPopSize,
               currentLogPopSize;

        // Get interval times
        prepareTimes();

        // Get likelihood for each segment
        logP = 0.0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {

            if (intervals.getIntervalType(i) == COALESCENT) {
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

        }

        return logP;
    }


    /**
     * Copied and simplified from BayesianSkyline.java
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
     * TODO check this!
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


    /************************************/
    /* Methods for getting change times */
    /************************************/

    public int getDimension() {
        return popSizes.getDimension();
    }


    /**
     * TODO: Fix it so it maps to the right coalescent intervals
     *
     * @param i
     * @return
     */
    public double getChangeTime(int i) {
        prepareTimes();

        int idx = Math.min(intervals.getIntervalCount()-1, cumulativeGroupSizes[i]);
        return coalescentTimes[idx];
    }


    /****************************/
    /* Calculation Node methods */
    /****************************/

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }


}
