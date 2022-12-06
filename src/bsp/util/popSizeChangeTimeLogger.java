package bsp.util;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;
import bsp.distributions.BSP;

import java.io.PrintStream;

public class popSizeChangeTimeLogger extends CalculationNode implements Loggable, Function {

    final public Input<BSP> skylineInput =
            new Input<>("skyline", "Skyline to log change times for", Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() {
        // nothing to do
    }

    @Override
    public void init(PrintStream out) {
        final BSP skyline = skylineInput.get();
        final int valueCount = skyline.getPopSizeDimension();

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
        final BSP skyline = skylineInput.get();

        final int values = skyline.getPopSizeDimension();
        for (int value = 0; value < values; value++) {
            out.print(skyline.getPopSizeChangeTime(value) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

    @Override
    public int getDimension() {
        final BSP skyline = skylineInput.get();
        return skyline.getDimension();
    }

    @Override
    public double getArrayValue() {
        final BSP skyline = skylineInput.get();
        return skyline.getPopSizeChangeTime(0);
    }

    @Override
    public double getArrayValue(int dim) {
        final BSP skyline = skylineInput.get();
        return skyline.getPopSizeChangeTime(dim);
    }



}
