package bsp.util;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;
import bsp.distributions.BESP;

import java.io.PrintStream;

public class samplingIntensityChangeTimeLogger extends CalculationNode implements Loggable, Function {

    final public Input<BESP> skylineInput =
            new Input<>("skyline", "Skyline to log change times for", Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() {
        // nothing to do
    }

    @Override
    public void init(PrintStream out) {
        final BESP skyline = skylineInput.get();
        final int valueCount = skyline.getSamplingIntensityDimension();

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
        final BESP skyline = skylineInput.get();

        final int values = skyline.getSamplingIntensityDimension();
        for (int value = 0; value < values; value++) {
            out.print(skyline.getSamplingIntensityChangeTime(value) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

    @Override
    public int getDimension() {
        final BESP skyline = skylineInput.get();
        return skyline.getDimension();
    }

    @Override
    public double getArrayValue() {
        final BESP skyline = skylineInput.get();
        return skyline.getSamplingIntensityChangeTime(0);
    }

    @Override
    public double getArrayValue(int dim) {
        final BESP skyline = skylineInput.get();
        return skyline.getSamplingIntensityChangeTime(dim);
    }

}
