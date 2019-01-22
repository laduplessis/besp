package bsp.distributions;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.coalescent.IntervalType;
import beast.math.Binomial;

import java.util.Arrays;

import static beast.evolution.tree.coalescent.IntervalType.COALESCENT;
import static beast.evolution.tree.coalescent.IntervalType.SAMPLE;


/**
 * Bayesian Skyline Plot implementation (with sampling intensity epochs/groups).
 *
 * @author Louis du Plessis
 * @date 2019/01/21
 */
public class PrefBSP extends BSP {

    final public Input<RealParameter> samplingIntensityInput =
            new Input<>("samplingIntensity","Constant sampling intensity");

    final public Input<IntegerParameter> samplingIntensityGroupSizeInput =
            new Input<>("samplingIntensityGroupSizes", "The number of events in each sampling intensity group in the skyline (use robust design if not provided)");


    protected RealParameter samplingIntensity;
    protected IntegerParameter samplingIntensityGroupSizes;

    protected int []    cumulativeSamplingIntensityGroupSizes;
    protected double [] samplingTimes, samplingIntensityGroupTimes;

    @Override
    public void initAndValidate() {

        int nrEvents, nrSamples,
            popGroups, samplingGroups;

        //////////////////////////
        // Get skyline parameters
        popSizes  = popSizeInput.get();
        popGroups = popSizes.getDimension();

        samplingIntensity = samplingIntensityInput.get();
        samplingGroups    = samplingIntensity.getDimension();

        ///////////////////////
        // Get tree intervals
        if (treeInput.get() != null) {
            throw new IllegalArgumentException("Only tree intervals (not tree) should be specified!");
        } else {
            intervals = treeIntervalsInput.get();
        }
        nrEvents  = intervals.getIntervalCount();
        nrSamples = intervals.getSampleCount()+1;

        ////////////////////////////
        // Get minimum group time
        minGroupTime = minGroupTimeInput.get();


        ////////////////////
        // Get group sizes
        // If group sizes are not specified use robust design (equal group sizes)
        if (popSizeGroupSizeInput.get() != null) {
            popSizeGroupSizes = popSizeGroupSizeInput.get();
        } else {
            popSizeGroupSizes = getRobustpopSizeGroupSizes(nrEvents, popGroups, 1, Integer.MAX_VALUE);
        }

        if (samplingIntensityGroupSizeInput.get() != null) {
            samplingIntensityGroupSizes = samplingIntensityGroupSizeInput.get();
        } else {
            samplingIntensityGroupSizes = getRobustpopSizeGroupSizes(nrSamples, samplingGroups, 1, Integer.MAX_VALUE);
        }

        // Group sizes does not equal the dimension of the skyline parameter
        if (popSizeGroupSizes.getDimension() != popGroups || samplingIntensityGroupSizes.getDimension() != samplingGroups) {
            throw new IllegalArgumentException("Number of groups should match the dimension of the skyline parameter "
                                             + "(effective population size or sampling intensity).");
        }

        // More groups than events
        if (popSizeGroupSizes.getDimension() > nrEvents || samplingIntensityGroupSizes.getDimension() > nrSamples) {
            throw new IllegalArgumentException("There are more groups than coalescent/sampling events in the tree.");
        }


        /////////////////////
        // Initialise arrays
        cumulativePopSizeGroupSizes           = new int[popGroups];
        cumulativeSamplingIntensityGroupSizes = new int[samplingGroups];
        //storedCumulativepopSizeGroupSizes = new int[nrGroups];
        popSizeGroupTimes           = new double[popGroups];
        samplingIntensityGroupTimes = new double[samplingGroups];
        updateArrays();

        // popSizeGroupSizes needs to add up to coalescent + sampling events
        if (cumulativePopSizeGroupSizes[popGroups - 1] != nrEvents) {
            Log.warning.println("WARNING: The sum of the initial effective population group sizes does not match the number of coalescent "
                              + "and sampling events in the tree. Initializing to equal group sizes (robust design)");

            popSizeGroupSizes.assignFromWithoutID(getRobustpopSizeGroupSizes(nrEvents, popGroups,
                                                  popSizeGroupSizes.getLower(), popSizeGroupSizes.getUpper()));

            // Recalculate cumulative group sizes, because group sizes have been changed
            updateArrays();
        }

        // samplingIntensityGroupSizes needs to add up to sampling events
        if (cumulativeSamplingIntensityGroupSizes[samplingGroups - 1] != nrSamples) {
            Log.warning.println("WARNING: The sum of the initial sampling intensity group sizes does not match the number of "
                              + "sampling events in the tree. Initializing to equal group sizes (robust design)");

            samplingIntensityGroupSizes.assignFromWithoutID(getRobustpopSizeGroupSizes(nrSamples, samplingGroups,
                    samplingIntensityGroupSizes.getLower(), samplingIntensityGroupSizes.getUpper()));

            // Recalculate cumulative group sizes, because group sizes have been changed
            updateArrays();
        }

    }



    /**
     * Updates the arrays used in likelihood calculation and other methods
     */
    protected void updateArrays() {

        // Get popsize cumulative group sizes and times (and extract sampling times)
        cumulativePopSizeGroupSizes[0] = popSizeGroupSizes.getValue(0);
        popSizeGroupTimes[0]           = intervals.getIntervalTime(cumulativePopSizeGroupSizes[0]-1);
        for (int i = 1; i < cumulativePopSizeGroupSizes.length; i++) {
            cumulativePopSizeGroupSizes[i] = cumulativePopSizeGroupSizes[i-1] + popSizeGroupSizes.getValue(i);
            popSizeGroupTimes[i]           = intervals.getIntervalTime(cumulativePopSizeGroupSizes[i]-1);
        }

        // Get sampling times
        samplingTimes = new double [intervals.getSampleCount()+1];
        int j = 0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            if (intervals.getIntervalType(i) == SAMPLE) {
                samplingTimes[j] = intervals.getIntervalTime(i);
                j++;
            }
        }

        // Get sampling intensity cumulative group sizes and times
        cumulativeSamplingIntensityGroupSizes[0] = samplingIntensityGroupSizes.getValue(0);
        samplingIntensityGroupTimes[0]           = samplingTimes[cumulativeSamplingIntensityGroupSizes[0]-1];
        for (int i = 1; i < cumulativeSamplingIntensityGroupSizes.length; i++) {
            cumulativeSamplingIntensityGroupSizes[i] = cumulativeSamplingIntensityGroupSizes[i-1] + samplingIntensityGroupSizes.getValue(i);
            samplingIntensityGroupTimes[i]           = samplingTimes[cumulativeSamplingIntensityGroupSizes[i]-1];
        }

        arraysUpdated = true;
    }



    @Override
    public double calculateLogP() {

        int    popSizeGroup           = 0,
               samplingIntensityGroup = 0,
               sampleIndex = 0;
        double width,
               currentTime = 0.0,
               currentPopSize,
               currentSamplingIntensity;

        // Update arrays
        if (!arraysUpdated) {
            updateArrays();

            //if (!(checkPopSizeGroupTimes() && checkPopSizeGroupSizes())) {
            //    return Double.NEGATIVE_INFINITY;
            //}
        }

        //System.out.println(Arrays.toString(samplingTimes));
        //System.out.println(Arrays.toString(cumulativeSamplingIntensityGroupSizes));

        // Get likelihood for each segment
        logP = 0.0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {

            // Next sampling intensity group
            if (intervals.getIntervalType(i) == SAMPLE) {
                sampleIndex++;
                if (sampleIndex > cumulativeSamplingIntensityGroupSizes[samplingIntensityGroup])
                    samplingIntensityGroup++;
            }

            // Next population size group
            if (i >= cumulativePopSizeGroupSizes[popSizeGroup]) {
                popSizeGroup++;
            }

            width = intervals.getInterval(i);
            currentTime += width;

            //System.out.println(samplingIntensityGroup);

            currentPopSize           = popSizes.getArrayValue(popSizeGroup);
            currentSamplingIntensity = currentTime <= samplingTimes[samplingTimes.length-1] ? samplingIntensity.getValue(samplingIntensityGroup) : 0.0;

            logP += calculateIntervalLikelihood(currentPopSize, currentSamplingIntensity, width, intervals.getLineageCount(i), intervals.getIntervalType(i));
        }

        return logP;
    }

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

        //System.out.println(lk + ":\t" + type+"\t"+lineageCount);
        //System.out.printf("%10.5f |%10.5f |%10d |%12s |%10s |%18s\n",
        //        lk, width, lineageCount, type, popSize, beta);


        return lk;
    }

}
