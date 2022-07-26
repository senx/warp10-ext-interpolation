//
// Copyright 2021 - 2022 SenX
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

package io.warp10;

import io.warp10.ext.interpolation.InterpolationWarpScriptExtension;
import io.warp10.script.MemoryWarpScriptStack;
import io.warp10.script.WarpScriptLib;
import io.warp10.script.WarpScriptStackFunction;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.junit.Assert;
import org.junit.Test;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;

public class CommonsMathWarpScriptExtensionTest {

  /**
   * Test case for the microsphere interpolator
   */
  @Test
  public void test() throws Exception {
    StringBuilder props = new StringBuilder();

    props.append("warp.timeunits=us\n");
    WarpConfig.safeSetProperties(new StringReader(props.toString()));
    WarpScriptLib.register(new InterpolationWarpScriptExtension());

    MemoryWarpScriptStack stack = new MemoryWarpScriptStack(null, null);
    stack.maxLimits();

    MultivariateFunction f = new MultivariateFunction() {
      @Override
      public double value(double[] x) {
        if (x.length != 2) {
          throw new IllegalArgumentException();
        }
        return 2 * x[0] - 3 * x[1] + 5;
      }
    };

    final double min = -1;
    final double max = 1;
    final double range = max - min;
    final int res = 5;
    final int n = res * res; // Number of sample points.
    final int dim = 2;

    List<List<Double>> x = new ArrayList<>(n);
    for (int i = 0; i < n; i++) {
      x.add(i, new ArrayList<Double>(dim));
    }
    List<Double> y = new ArrayList<>(n);

    int index = 0;
    for (int i = 0; i < res; i++) {
      final double x1Val = toCoordinate(min, range, res, i);
      for (int j = 0; j < res; j++) {
        final double x2Val = toCoordinate(min, range, res, j);

        x.get(index).add(x1Val);
        x.get(index).add(x2Val);

        double[] point = new double[dim];
        for (int k = 0; k < dim; k++) {
          point[k] = x.get(index).get(k);
        }
        y.add(f.value(point));

        ++index;
      }
    }

    stack.push(x);
    stack.push(y);
    stack.execMulti("MICROSPHEREFIT");
    Object func = stack.pop();

    double[] c = new double[dim];
    double expected;
    double result;

    final int sampleIndex = 2;
    c[0] = x.get(sampleIndex).get(0);
    c[1] = x.get(sampleIndex).get(1);
    expected = f.value(c);
    List cList = new ArrayList<Double>(2);
    cList.add(c[0]);
    cList.add(c[1]);
    stack.push(cList);
    ((WarpScriptStackFunction) func).apply(stack);
    result = (Double) stack.pop();
    Assert.assertEquals("on sample point (exact)", expected, result, Math.ulp(1d));

    // Interpolation.
    c[0] = 0.654321;
    c[1] = -0.345678;
    expected = f.value(c);
    cList = new ArrayList<Double>(2);
    cList.add(c[0]);
    cList.add(c[1]);
    stack.push(cList);
    ((WarpScriptStackFunction) func).apply(stack);
    result = (Double) stack.pop();
    Assert.assertEquals("interpolation (exact)", expected, result, 1e-1);

    // Extrapolation.
    c[0] = 0 - 1e-2;
    c[1] = 1 + 1e-2;
    expected = f.value(c);
    cList = new ArrayList<Double>(2);
    cList.add(c[0]);
    cList.add(c[1]);
    stack.push(cList);
    ((WarpScriptStackFunction) func).apply(stack);
    result = (Double) stack.pop();
    Assert.assertFalse(Double.isNaN(result));
    Assert.assertEquals("extrapolation (exact)", expected, result, 1e-1);

    // Far away.
    c[0] = 20;
    c[1] = -30;
    cList = new ArrayList<Double>(2);
    cList.add(c[0]);
    cList.add(c[1]);
    stack.push(cList);
    ((WarpScriptStackFunction) func).apply(stack);
    result = (Double) stack.pop();
    Assert.assertTrue(result + " should be NaN", Double.isNaN(result));
  }

  /**
   * @param min Minimum of the coordinate range.
   * @param range Extent of the coordinate interval.
   * @param res Number of pixels.
   * @param pixel Pixel index.
   */
  private static double toCoordinate(double min, double range, int res, int pixel) {
    return pixel * range / (res - 1) + min;
  }

}