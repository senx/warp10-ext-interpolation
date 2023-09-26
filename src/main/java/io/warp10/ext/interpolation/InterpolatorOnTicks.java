//
//   Copyright 2023  SenX S.A.S.
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
//

package io.warp10.ext.interpolation;

import io.warp10.continuum.gts.GTSHelper;
import io.warp10.continuum.gts.GeoTimeSerie;
import io.warp10.script.NamedWarpScriptFunction;
import io.warp10.script.WarpScriptException;
import io.warp10.script.WarpScriptMapperFunction;
import io.warp10.script.WarpScriptStack;
import io.warp10.script.WarpScriptStackFunction;
import io.warp10.script.functions.SNAPSHOT;
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

/**
 * Function that implements 1D interpolator on ticks
 */
public class InterpolatorOnTicks extends NamedWarpScriptFunction implements WarpScriptStackFunction {

  public enum TYPE {
    LINEAR, SPLINE, AKIMA
  }

  private TYPE type;

  private static class MAPPER_INTERPOLATOR extends NamedWarpScriptFunction implements WarpScriptMapperFunction {

    private final PolynomialSplineFunction func;
    private final String generatedFrom;
    private GeoTimeSerie gts;

    private MAPPER_INTERPOLATOR(PolynomialSplineFunction function, String interpolatorName, GeoTimeSerie gts) {
      super(interpolatorName);
      func = function;
      this.gts = gts;
      generatedFrom = interpolatorName;
    }

    private double value(double x) {
      if (!func.isValidPoint(x)) {
        return Double.NaN;
      } else {
        return func.value(x);
      }
    }

    @Override
    public Object apply(Object[] args) throws WarpScriptException {
      long tick = (long) args[0];
      long[] ticks = (long[]) args[3];
      long[] locations = (long[]) args[4];
      long[] elevations = (long[]) args[5];
      Object[] values = (Object[]) args[6];

      if (1 < ticks.length) {
        throw new WarpScriptException(getName() + " expects at most 1 tick but got " + ticks.length);
      }

      //
      // We return the same value if a value is present,
      // else, we interpolate.
      //

      if (1 == values.length) {
        double value = ((Number) values[0]).doubleValue();

        if (!Double.isNaN(value)) {
          return new Object[] {tick, locations[0], elevations[0], value};

        } else {
          return new Object[] {tick, locations[0], elevations[0], value(tick)};
        }
      }

      return new Object[] {tick, locations[0], elevations[0], value(tick)};
    }

    @Override
    public String toString() {
      StringBuilder sb = new StringBuilder();
      try {
        SNAPSHOT.addElement(sb, gts);
        sb.append(" ");
      } catch (Exception e) {
        throw new RuntimeException("Error building argument snapshot", e);
      }

      sb.append(generatedFrom);

      return sb.toString();
    }
  }

  public InterpolatorOnTicks(String name, TYPE type) {
    super(name);
    this.type=type;
  }

  @Override
  public Object apply(WarpScriptStack stack) throws WarpScriptException {

    double xval[];
    double fval[];

    Object o = stack.pop();

    if (!(o instanceof GeoTimeSerie)) {
      throw new WarpScriptException(getName() + " expects a GTS");
    }
    GeoTimeSerie gts = (GeoTimeSerie) o;
    
    // ticks need to be sorted in order to create the interpolating function
    GTSHelper.sort(gts);

    // fill x and f
    int size = gts.size();
    if (getName().equals("AKIMASPLINEFIT") && size < 5) {
      throw new WarpScriptException(getName() + " expects at least 5 interpolation points");
    }
    xval = new double[size];
    fval = new double[size];
    for (int i = 0; i < size; i++) {
      xval[i] = ((Number) GTSHelper.tickAtIndex(gts, i)).doubleValue();
      fval[i] = ((Number) GTSHelper.valueAtIndex(gts, i)).doubleValue();
    }

    PolynomialSplineFunction function = null;
    if (type== TYPE.LINEAR) {
      function = (new LinearInterpolator()).interpolate(xval, fval);
    } else if (type== TYPE.SPLINE) {
      function = (new SplineInterpolator()).interpolate(xval, fval);
    } else if (type== TYPE.AKIMA) {
      function = (new AkimaSplineInterpolator()).interpolate(xval, fval);
    }

    // clone the inputs for snapshot. 
    MAPPER_INTERPOLATOR warpscriptFunction = new MAPPER_INTERPOLATOR(function, getName(), gts.clone());
    stack.push(warpscriptFunction);

    return stack;
  }
}
