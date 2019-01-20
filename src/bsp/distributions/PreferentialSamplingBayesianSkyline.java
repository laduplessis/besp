package bsp.distributions;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.coalescent.IntervalType;
import beast.math.Binomial;

import static beast.evolution.tree.coalescent.IntervalType.COALESCENT;


/**
 * Basic preferential sampling Bayesian skyline implementation
 *
 * This implementation only allows linear preferential sampling
 * In addition, the sampling intensity is assumed to be constant through time
 *
 * The major difference to the Bayesian skyline plot is that groups are delineated not only at coalescent events, but
 * also at sampling events (which means that it mostly needs to be rewritten)
 *
 *
 */
public class PreferentialSamplingBayesianSkyline extends GeneralBayesianSkyline {

    final public Input<RealParameter> samplingIntensityInput =
            new Input<>("samplingIntensity","Constant sampling intensity");


    protected RealParameter samplingIntensity;
    protected double oldestSample;

    @Override
    public void initAndValidate() {

        int nrEvents, nrGroups;

        //////////////////////////
        // Get skyline parameter
        if (popSizeInput.get() == null) {
            useLog = true;
            logPopSizes = logPopSizeInput.get();
            nrGroups    = logPopSizes.getDimension();
            throw new IllegalArgumentException("Log population size parameterization not tested yet!");
        } else {
            useLog = false;
            popSizes = popSizeInput.get();
            nrGroups = popSizes.getDimension();
        }

        //////////////////////////
        // Get sampling intensity
        samplingIntensity = samplingIntensityInput.get();

        ///////////////////////
        // Get tree intervals
        if (treeInput.get() != null) {
            // Doesn't work... tree intervals are not updated with tree
            // intervals.initByName("tree", treeInput.get());
            throw new IllegalArgumentException("Only tree intervals (not tree) should be specified!");
        } else {
            intervals = treeIntervalsInput.get();
        }
        nrEvents = intervals.getIntervalCount();


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

        cumulativeGroupSizes = new int[nrGroups];
        groupTimes = new double[nrGroups];
        updateArrays();

        // GroupSizes needs to add up to coalescent events
        if (cumulativeGroupSizes[nrGroups - 1] != nrEvents) {
            Log.warning.println("WARNING: The sum of the initial group sizes does not match the number of coalescent " +
                    "and sampling events in the tree. Initializing to equal group sizes (robust design)");

            groupSizes.assignFromWithoutID(getRobustGroupSizes(nrEvents, nrGroups));

            // Recalculate cumulative group sizes, because group sizes have been changed
            updateArrays();
        }
    }


    /**
     * Updates the arrays used in likelihood calculation and other methods
     */
    @Override
    protected void updateArrays() {

        // Get coalescent times (not necessary for likelihood computation)
        // coalescentTimes = intervals.getCoalescentTimes(coalescentTimes);

        // Get cumulative group sizes and times
        cumulativeGroupSizes[0] = groupSizes.getValue(0);
        groupTimes[0]           = intervals.getIntervalTime(cumulativeGroupSizes[0]-1);
        for (int i = 1; i < cumulativeGroupSizes.length; i++) {
            cumulativeGroupSizes[i] = cumulativeGroupSizes[i-1] + groupSizes.getValue(i);
            groupTimes[i]           = intervals.getIntervalTime(cumulativeGroupSizes[i]-1);
        }

        // Get oldest sampling time (per epoch)
        // TODO can be improved
        double currentTime = 0.0;
        oldestSample = 0.0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            currentTime += intervals.getInterval(i);

            if (intervals.getIntervalType(i) == IntervalType.SAMPLE) {
                oldestSample = currentTime;
            }
        }

        arraysUpdated = true;
    }


    @Override
    public double calculateLogP() {

        int    groupIndex = 0;
        double width,
               currentTime = 0.0,
               currentPopSize,
               currentsamplingIntensity;

        // Update arrays
        if (!arraysUpdated) {
            updateArrays();
        }

        // Get likelihood for each segment
        logP = 0.0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {

            if (i >= cumulativeGroupSizes[groupIndex])
                groupIndex++;

            width = intervals.getInterval(i);
            currentTime += width;

            currentPopSize = popSizes.getArrayValue(groupIndex);
            currentsamplingIntensity = currentTime <= oldestSample ? samplingIntensity.getValue() : 0.0;


            //System.out.println(i+"\t"+coalIndex+"\t"+groupIndex+"\t"+currentTime+"\t"+intervals.getInterval(i)+"\t"+currentPopSize+"\t"+currentPopSize2);
            //System.out.println(i+"\t"+coalIndex+"\t"+groupIndex+"\t"+currentTime+"\t"+intervals.getInterval(i)+"\t"+intervals.getCoalescentEvents(i));

            // Calculate likelihood for interval
            logP += calculateIntervalLikelihood(currentPopSize, currentsamplingIntensity, width, intervals.getLineageCount(i), intervals.getIntervalType(i));



            //System.out.println(logP + "\t" + currentTime);

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
        //System.out.printf("%10.5f %10s |%10d%12s |%10s%18s\n",
        //        lk, width, lineageCount, type, popSize, beta);


        return lk;
    }


}
