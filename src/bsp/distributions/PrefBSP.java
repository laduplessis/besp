package bsp.distributions;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.coalescent.IntervalType;
import beast.math.Binomial;
import org.omg.SendingContext.RunTime;

import java.util.Arrays;

import static beast.evolution.tree.coalescent.IntervalType.SAMPLE;


/***************************************************************
 * Bayesian Epoch Sampling Skyline Plot implementation (BESP). *
 * *************************************************************
 *
 * @author Louis du Plessis
 * @date 2019/01/21
 */
@Description("BESP: A likelihood function for epoch sampling skyline plot coalescent.")
@Citation(value="Parag, K.V., du Plessis, L., Pybus, O.G. (2019).\n"+
        "  Jointly inferring the dynamics of population size and sampling intensity from molecular sequences.\n",
        year = 2019, firstAuthorSurname = "Parag", DOI="10.1101/686378")

public class PrefBSP extends BSP {

    final public Input<RealParameter> samplingIntensityInput =
            new Input<>("samplingIntensity","Sampling intensity (for each epoch)");

    final public Input<IntegerParameter> samplingIntensityGroupSizeInput =
            new Input<>("samplingIntensityGroupSizes", "The number of sampling events in each sampling intensity epoch (use robust design if not provided)");

    final public Input<RealParameter> samplingEpochTimesInput =
            new Input<>("samplingEpochTimes","Times when the sampling intensity change (distance from most recent tip)");

    protected RealParameter samplingIntensity;
    protected IntegerParameter samplingIntensityGroupSizes;

    protected int []    cumulativeSamplingIntensityGroupSizes;
    protected double [] samplingTimes, intervalTimes,
                        samplingIntensityGroupTimes;

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
        minWidth = minWidthInput.get();


        ////////////////////
        // Get group sizes
        // If epoch times are specified use those to get the groups
        // If group sizes are not specified use robust design (equal group sizes)

        // popSize
        if (popSizeEpochTimesInput.get() != null) {

            if (popSizeGroupSizeInput.get() != null) {
                throw new IllegalArgumentException("Only one of popSizeEpochTimes and popSizeGroupSizes should be specified.");
            }

            RealParameter popSizeEpochTimes = popSizeEpochTimesInput.get();
            double [] eventTimes = new double [intervals.getIntervalCount()];
            for (int i = 0; i < eventTimes.length; i++) {
                eventTimes[i] = intervals.getIntervalTime(i);
            }
            popSizeGroupSizes = epochsToGroups(eventTimes, popSizeEpochTimes.getValues(), 1, Integer.MAX_VALUE);

        } else {

            if (popSizeGroupSizeInput.get() != null) {
                popSizeGroupSizes = popSizeGroupSizeInput.get();
            } else {
                popSizeGroupSizes = getRobustPopSizeGroupSizes(nrEvents, popGroups, 1, Integer.MAX_VALUE);
            }

        }

        // samplingIntensity
        if (samplingEpochTimesInput.get() != null) {

            if (samplingIntensityGroupSizeInput.get() != null) {
                throw new IllegalArgumentException("Only one of samplingEpochTimes and samplingIntensityGroupSizes should be specified.");
            }

            RealParameter samplingEpochTimes = samplingEpochTimesInput.get();

            intervalTimes = new double [intervals.getIntervalCount()];
            samplingTimes = new double [intervals.getSampleCount()+1];
            updateIntervalTimes(intervalTimes, samplingTimes);
            samplingIntensityGroupSizes = epochsToGroups(samplingTimes, samplingEpochTimes.getValues(), 1, Integer.MAX_VALUE);

        } else {

            if (samplingIntensityGroupSizeInput.get() != null) {
                samplingIntensityGroupSizes = samplingIntensityGroupSizeInput.get();
            } else {
                samplingIntensityGroupSizes = getRobustPopSizeGroupSizes(nrSamples, samplingGroups, 1, Integer.MAX_VALUE);
            }

        }

        // Group sizes do not equal the dimension of the skyline parameter
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
        intervalTimes = new double [intervals.getIntervalCount()];
        samplingTimes = new double [intervals.getSampleCount()+1];
        //samplingTimes = getSamplingTimes(getIntervalTimes(null));
        updateArrays();

        // popSizeGroupSizes needs to add up to coalescent + sampling events
        if (cumulativePopSizeGroupSizes[popGroups - 1] != nrEvents) {
            Log.warning.println("WARNING: The sum of the initial effective population group sizes does not match the number of coalescent "
                              + "and sampling events in the tree. Initializing to equal group sizes (robust design)");

            popSizeGroupSizes.assignFromWithoutID(getRobustPopSizeGroupSizes(nrEvents, popGroups,
                                                  popSizeGroupSizes.getLower(), popSizeGroupSizes.getUpper()));

            // Recalculate cumulative group sizes, because group sizes have been changed
            updateArrays();
        }

        // samplingIntensityGroupSizes needs to add up to sampling events
        if (cumulativeSamplingIntensityGroupSizes[samplingGroups - 1] != nrSamples) {
            Log.warning.println("WARNING: The sum of the initial sampling intensity group sizes does not match the number of "
                              + "sampling events in the tree. Initializing to equal group sizes (robust design)");

            samplingIntensityGroupSizes.assignFromWithoutID(getRobustPopSizeGroupSizes(nrSamples, samplingGroups,
                    samplingIntensityGroupSizes.getLower(), samplingIntensityGroupSizes.getUpper()));

            // Recalculate cumulative group sizes, because group sizes have been changed
            updateArrays();
        }

        // popSize group widths need to be longer than minWidth
        int i = 0;
        int numInitializationAttemps = numInitializationAttemptsInput.get();
        while (!checkGroupWidths(popSizeGroupTimes, minWidth)) {
            if (i == 0) {
                Log.warning.println("WARNING: Minimum effective population group width is shorter than minWidth ("+ minWidth + ")\n"
                                  + "Attempting to adjust...");
            } else
            if (i > numInitializationAttemps) {
                throw new IllegalArgumentException("Minimum effective population group width is still shorter than minWidth (" + minWidth + ") "
                                                 + "after "+numInitializationAttemps+" attempts at redistributing group sizes.\n"
                                                 + "Try decreasing the number of groups or the minimum group width.");
                                               //+ "Current group times: " + Arrays.toString(popSizeGroupTimes));
            }
            redistributeGroups(popSizeGroupSizes, popSizeGroupTimes);
            updateArrays();
            i++;
        }

        // samplingIntensity group widths need to be longer than minWidth
        i = 0;
        while (!checkGroupWidths(samplingIntensityGroupTimes, minWidth)) {
            if (i == 0) {
                Log.warning.println("WARNING: Minimum  sampling intensity group width is shorter than minWidth ("+ minWidth + ")\n"
                                  + "Attempting to adjust...");
            } else
            if (i > numInitializationAttemps) {
                throw new IllegalArgumentException("Minimum sampling intensity group width is still shorter than minWidth (" + minWidth + ") "
                                                 + "after "+numInitializationAttemps+" attempts at redistributing group sizes.\n"
                                                 + "Try decreasing the number of groups or the minimum group width.");
                                               //+ "Current group times: " + Arrays.toString(popSizeGroupTimes));
            }
            redistributeGroups(samplingIntensityGroupSizes, samplingIntensityGroupTimes);
            updateArrays();
            i++;
        }

        System.out.println(this.toString());

    }



    protected void updateIntervalTimes(double [] intervalTimes, double [] samplingTimes) {

        double time = 0.0;
        int j = 0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            time += intervals.getInterval(i);
            intervalTimes[i] = time;
            if (intervals.getIntervalType(i) == SAMPLE) {
                samplingTimes[j] = time;
                j++;
            }
        }
    }


    /**
     * Updates the arrays used in likelihood calculation and other methods
     */
    protected void updateArrays() {

        //System.out.println("Updating arrays...");
        updateIntervalTimes(intervalTimes, samplingTimes);

        // Get popsize cumulative group sizes and times (and extract sampling times)
        cumulativePopSizeGroupSizes[0] = popSizeGroupSizes.getValue(0);
        popSizeGroupTimes[0]           = intervalTimes[cumulativePopSizeGroupSizes[0]-1];
        for (int i = 1; i < cumulativePopSizeGroupSizes.length; i++) {
            cumulativePopSizeGroupSizes[i] = cumulativePopSizeGroupSizes[i-1] + popSizeGroupSizes.getValue(i);
            popSizeGroupTimes[i]           = intervalTimes[cumulativePopSizeGroupSizes[i]-1];
        }

        //System.out.println((cumulativePopSizeGroupSizes[cumulativePopSizeGroupSizes.length-1]-1)+"\t"+(intervals.getIntervalCount()-1));

        //System.out.println("grouptime: "+popSizeGroupTimes[popSizeGroupTimes.length-1]+"\ttree intervals: "+
        //        intervals.getIntervalTime(intervals.getIntervalCount()-1) +"\t"+
        //        (popSizeGroupTimes[popSizeGroupTimes.length-1] == intervals.getIntervalTime(intervals.getIntervalCount()-1)));

        // Get sampling times
        //samplingTimes = getSamplingTimes(intervalTimes);

        // Get sampling intensity cumulative group sizes and times
        cumulativeSamplingIntensityGroupSizes[0] = samplingIntensityGroupSizes.getValue(0);
        samplingIntensityGroupTimes[0]           = samplingTimes[cumulativeSamplingIntensityGroupSizes[0]-1];
        for (int i = 1; i < cumulativeSamplingIntensityGroupSizes.length; i++) {
            cumulativeSamplingIntensityGroupSizes[i] = cumulativeSamplingIntensityGroupSizes[i-1] + samplingIntensityGroupSizes.getValue(i);
            samplingIntensityGroupTimes[i]           = samplingTimes[cumulativeSamplingIntensityGroupSizes[i]-1];
        }




        /*
        // Get popsize cumulative group sizes and times (and extract sampling times)
        cumulativePopSizeGroupSizes[0] = popSizeGroupSizes.getValue(0);
        popSizeGroupTimes[0]           = intervals.getIntervalTime(cumulativePopSizeGroupSizes[0]-1);
        for (int i = 1; i < cumulativePopSizeGroupSizes.length; i++) {
            cumulativePopSizeGroupSizes[i] = cumulativePopSizeGroupSizes[i-1] + popSizeGroupSizes.getValue(i);
            popSizeGroupTimes[i]           = intervals.getIntervalTime(cumulativePopSizeGroupSizes[i]-1);
        }

        //System.out.println((cumulativePopSizeGroupSizes[cumulativePopSizeGroupSizes.length-1]-1)+"\t"+(intervals.getIntervalCount()-1));

        //System.out.println("grouptime: "+popSizeGroupTimes[popSizeGroupTimes.length-1]+"\ttree intervals: "+
        //        intervals.getIntervalTime(intervals.getIntervalCount()-1) +"\t"+
        //        (popSizeGroupTimes[popSizeGroupTimes.length-1] == intervals.getIntervalTime(intervals.getIntervalCount()-1)));

        // Get sampling times
        samplingTimes = getSamplingTimes();

        // Get sampling intensity cumulative group sizes and times
        cumulativeSamplingIntensityGroupSizes[0] = samplingIntensityGroupSizes.getValue(0);
        samplingIntensityGroupTimes[0]           = samplingTimes[cumulativeSamplingIntensityGroupSizes[0]-1];
        for (int i = 1; i < cumulativeSamplingIntensityGroupSizes.length; i++) {
            cumulativeSamplingIntensityGroupSizes[i] = cumulativeSamplingIntensityGroupSizes[i-1] + samplingIntensityGroupSizes.getValue(i);
            samplingIntensityGroupTimes[i]           = samplingTimes[cumulativeSamplingIntensityGroupSizes[i]-1];
        }
         */

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

        //System.out.println("Calculating likelihood...");

        // Update arrays
        if (!arraysUpdated) {
            updateArrays();

            if (!(checkGroupWidths(popSizeGroupTimes, minWidth) && checkGroupWidths(samplingIntensityGroupTimes, minWidth))) {
                return Double.NEGATIVE_INFINITY;
            }
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

    @Override
    public String toString() {

        double start  = 0.0;
        String outstr = super.toString();

        outstr += String.format("%10s  %10s | %7s  %7s  %7s | %10s\n"+
                        "------------------------------------------------------------------------------\n",
                "group", "size", "start", "end", "width", "samplingIntensity");

        for (int i = 0; i < samplingIntensityGroupSizes.getDimension(); i++) {
            outstr += String.format("%10s  %10s | %4.5f  %4.5f  %4.5f | %10s\n",
                    i+1, samplingIntensityGroupSizes.getValue(i), start, samplingIntensityGroupTimes[i],
                    samplingIntensityGroupTimes[i]-start, samplingIntensity.getValue(i));
            start = samplingIntensityGroupTimes[i];
        }

        return outstr+"\n";

    }


    /**********************************/
    /* Methods for sampling intensity */
    /**********************************/

    public int getSamplingIntensityDimension() { return samplingIntensity.getDimension(); }

    /**
     * Return the i'th change time of the skyline parameter (for logging and skyline reconstruction)
     *
     * @param i
     * @return
     */
    public double getSamplingIntensityChangeTime(int i) {
        if (!arraysUpdated) {
            updateArrays();
        }

        return samplingIntensityGroupTimes[i];
    }


    /**
     * TODO: Check boundary conditions
     *
     * @param t
     * @return
     */
    public double getSamplingIntensity(double t) {

        int groupIndex = Arrays.binarySearch(samplingIntensityGroupTimes, t);
        if (groupIndex < 0) {
            groupIndex = -groupIndex - 1;
        } else {
            groupIndex++;
        }

        return popSizes.getValue(Math.min(groupIndex, samplingIntensity.getDimension()-1));
    }




}
