package bsp.util;

import beast.core.CalculationNode;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;
import bsp.distributions.GeneralBayesianSkyline;

import java.io.PrintStream;

public class BayesianSkylineChangeTimeLogger extends CalculationNode implements Loggable, Function {

    final public Input<GeneralBayesianSkyline> skylineInput =
            new Input<>("skyline", "Skyline to log change times for", Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() {
        // nothing to do
    }

    @Override
    public void init(PrintStream out) {
        final GeneralBayesianSkyline skyline = skylineInput.get();
        final int valueCount = skyline.getDimension();

        if (valueCount == 1) {
            out.print(this.getID()+"\t");
        } else {
            for (int value = 0; value < valueCount; value++) {
                out.print(this.getID()+ (value + 1) + "\t");
            }
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        final GeneralBayesianSkyline skyline = skylineInput.get();

        final int values = skyline.getDimension();
        for (int value = 0; value < values; value++) {
            out.print(skyline.getChangeTime(value) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

    @Override
    public int getDimension() {
        final GeneralBayesianSkyline skyline = skylineInput.get();
        return skyline.getDimension();
    }

    @Override
    public double getArrayValue() {
        final GeneralBayesianSkyline skyline = skylineInput.get();
        return skyline.getChangeTime(0);
    }

    @Override
    public double getArrayValue(int dim) {
        final GeneralBayesianSkyline skyline = skylineInput.get();
        return skyline.getChangeTime(dim);
    }



}
