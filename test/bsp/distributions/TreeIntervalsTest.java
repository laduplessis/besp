package bsp.distributions;

import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.IntervalType;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.util.TreeParser;
import junit.framework.TestCase;
import org.junit.Test;

import static beast.evolution.tree.coalescent.IntervalType.SAMPLE;

public class TreeIntervalsTest extends TestCase {


    @Test
    public void test1 () {

        // Tree with 6 samples, 2 sampling events at time 0 and at time 6, all other samples and coalescences at unique times
        // Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:19.0,(D4Thai78:5.0,D4Thai84:11.0):14.0):15.0);",false);

        // Tree with 6 samples, 2 sampling events at time 0 and at time 6 coalescences at time 23
        // Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:17.0,(D4Thai78:5.0,D4Thai84:11.0):12.0):17.0);",false);

        // Tree with 7 samples, polytomy at time 11
        Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:5.0,D4Thai78:5.0,D4Thai84:11.0,D4Test:2.0):29.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        System.out.println("Intervals: "+intervals.getIntervalCount());
        System.out.println("Coalescences: "+intervals.getSampleCount());

        double currentTime = 0.0;
        double width;
        String type;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            width = intervals.getInterval(i);
            type  = intervals.getIntervalType(i) == IntervalType.COALESCENT ? "coalescent" : "sampling";

            System.out.println(currentTime+"\t"+width+"\t"+intervals.getLineageCount(i)+"\t"+type+"\t"+intervals.getCoalescentEvents(i));
            currentTime += width;
        }

    }


    @Test
    public void testIntervalTimes () {

        // Tree with 6 samples, 2 sampling events at time 0 and at time 6, all other samples and coalescences at unique times
        //      Sampling times:  0, 0, 6, 6, 20, 28
        //      Branching times: 11, 23, 25, 30, 40
        Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:19.0,(D4Thai78:5.0,D4Thai84:11.0):14.0):15.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        double [] intervalTimes = new double[intervals.getIntervalCount()];

        double time = 0.0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            time += intervals.getInterval(i);
            intervalTimes[i] = time;
        }

        double [] expected = {0, 0, 6, 6, 11, 20, 23, 25, 28, 30, 40};
        for (int i = 0; i < intervalTimes.length; i++) {
            System.out.println(String.format("%5s\t%5s\t%10s", expected[i], intervalTimes[i], intervals.getIntervalType(i)));
            assertEquals(expected[i], intervalTimes[i]);
        }

    }

    @Test
    public void testSamplingTimes1 () {

        // Tree with 6 samples, 2 sampling events at time 0 and at time 6, all other samples and coalescences at unique times
        //      Sampling times:  0, 0, 6, 6, 20, 28
        //      Branching times: 11, 23, 25, 30, 40
        Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:19.0,(D4Thai78:5.0,D4Thai84:11.0):14.0):15.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        double [] values = new double [intervals.getSampleCount()+1];
        double [] expected = {0, 0, 6, 6, 20, 28};

        int j = 0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            if (intervals.getIntervalType(i) == SAMPLE) {
                values[j] = intervals.getIntervalTime(i);
                System.out.println(String.format("%5s\t%5s\t%10s", expected[j], values[j], intervals.getIntervalType(i)));
                j++;
            }
        }

        for (int i = 0; i < values.length; i++) {
            //System.out.println(String.format("%5s\t%5s\t%10s", expected[i], values[i], intervals.getIntervalType(i)));
            assertEquals(expected[i], values[i]);
        }

    }

    @Test
    public void testSamplingTimes2 () {

        // Tree with 6 samples, 2 sampling events at time 0 and at time 6, all other samples and coalescences at unique times
        //      Sampling times:  0, 2, 6, 6, 20, 28
        Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:19.0,(D4Thai78:5.0,D4Thai84:9.0):14.0):15.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        double [] values = new double [intervals.getSampleCount()+1];
        double [] expected = {0, 2, 6, 6, 20, 28};

        int j = 0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            if (intervals.getIntervalType(i) == SAMPLE) {
                values[j] = intervals.getIntervalTime(i);
                System.out.println(String.format("%5s\t%5s\t%10s", expected[j], values[j], intervals.getIntervalType(i)));
                j++;
            } else {
                System.out.println(String.format("%5s\t%5s\t%10s", "-", intervals.getIntervalTime(i), intervals.getIntervalType(i)));
            }

        }

        for (int i = 0; i < values.length; i++) {
            //System.out.println(String.format("%5s\t%5s\t%10s", expected[i], values[i], intervals.getIntervalType(i)));
            assertEquals(expected[i], values[i]);
        }
    }


    @Test
    public void testSamplingTimes3 () {

        //      Sampling times:  0, 0 6, 6, 9, 20, 28
        Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:5.0,D4Thai78:5.0,D4Thai84:11.0,D4Test:2.0):29.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        double [] values = new double [intervals.getSampleCount()+1];
        double [] expected = {0, 0, 6, 6, 9, 20, 28};

        int j = 0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            if (intervals.getIntervalType(i) == SAMPLE) {
                values[j] = intervals.getIntervalTime(i);
                System.out.println(String.format("%5s\t%5s\t%10s", expected[j], values[j], intervals.getIntervalType(i)));
                j++;
            } else {
                System.out.println(String.format("%5s\t%5s\t%10s", "-", intervals.getIntervalTime(i), intervals.getIntervalType(i)));
            }

        }

        for (int i = 0; i < values.length; i++) {
            //System.out.println(String.format("%5s\t%5s\t%10s", expected[i], values[i], intervals.getIntervalType(i)));
            assertEquals(expected[i], values[i]);
        }
    }

}
