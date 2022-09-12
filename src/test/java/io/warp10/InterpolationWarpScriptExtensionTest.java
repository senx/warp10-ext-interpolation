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
import io.warp10.ext.interpolation.MICROSPHEREFIT;
import io.warp10.ext.interpolation.TRICUBICFIT;
import io.warp10.script.MemoryWarpScriptStack;
import io.warp10.script.WarpScriptLib;
import io.warp10.script.WarpScriptStackFunction;
import io.warp10.script.ext.inventory.InventoryWarpScriptExtension;
import io.warp10.script.ext.token.TokenWarpScriptExtension;
import io.warp10.standalone.Warp;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.TrivariateFunction;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.StringReader;
import java.lang.reflect.Method;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class InterpolationWarpScriptExtensionTest {

  @Test
  public void tricubic_testPlane()  throws Exception {
    StringBuilder props = new StringBuilder();

    props.append("warp.timeunits=us\n");
    WarpConfig.safeSetProperties(new StringReader(props.toString()));
    InterpolationWarpScriptExtension ext = new InterpolationWarpScriptExtension();
    WarpScriptLib.register(ext);

    MemoryWarpScriptStack stack = new MemoryWarpScriptStack(null, null);
    stack.maxLimits();
    double[] xval = new double[] {3, 4, 5, 6.5};
    double[] yval = new double[] {-4, -3, -1, 2, 2.5};
    double[] zval = new double[] {-12, -8, -5.5, -3, 0, 2.5};

    List<Double> xlist = new ArrayList<>();
    for (int i = 0; i < xval.length; i++) {
      xlist.add(xval[i]);
    }
    List<Double> ylist = new ArrayList<>();
    for (int i = 0; i < yval.length; i++) {
      ylist.add(yval[i]);
    }
    List<Double> zlist = new ArrayList<>();
    for (int i = 0; i < zval.length; i++) {
      zlist.add(zval[i]);
    }

    // Function values
    TrivariateFunction f = new TrivariateFunction() {
      @Override
      public double value(double x, double y, double z) {
        return 2 * x - 3 * y - 4 * z + 5;
      }
    };

    double[][][] fval = new double[xval.length][yval.length][zval.length];
    List<List<List<Double>>> fvalTensor = new ArrayList<List<List<Double>>>();

    for (int i = 0; i < xval.length; i++) {
      fvalTensor.add(new ArrayList<>());

      for (int j = 0; j < yval.length; j++) {
        fvalTensor.get(i).add(new ArrayList<Double>());

        for (int k = 0; k < zval.length; k++) {
          // tensor
          List<Double> l = fvalTensor.get(i).get(j);
          l.add(f.value(xval[i], yval[j], zval[k]));

          // 3d arrays
          fval[i][j][k] = f.value(xval[i], yval[j], zval[k]);
        }
      }
    }

    stack.push(xlist);
    stack.push(ylist);
    stack.push(zlist);
    stack.push(fvalTensor);

    ((WarpScriptStackFunction) ext.getFunctions().get(TRICUBICFIT.class.getSimpleName())).apply(stack);
    WarpScriptStackFunction function = (WarpScriptStackFunction) stack.pop();

    List<Number> input = new ArrayList<>();
    input.add(4L);
    input.add(-3L);
    input.add(0L);
    stack.push(input);
    function.apply(stack);
    stack.push(f.value(input.get(0).doubleValue(), input.get(1).doubleValue(), input.get(2).doubleValue()));
    stack.execMulti("== ASSERT");

    input = new ArrayList<>();
    input.add(4.5);
    input.add(-1.5);
    input.add(-4.25);
    stack.push(input);
    function.apply(stack);
    stack.push(f.value(input.get(0).doubleValue(), input.get(1).doubleValue(), input.get(2).doubleValue()));
    stack.execMulti("== ASSERT");

    System.out.println(stack.dump(100));
  }

  @Test
  public void tricubic_testPrPFailCase()  throws Exception {
    StringBuilder props = new StringBuilder();

    props.append("warp.timeunits=us\n");
    WarpConfig.safeSetProperties(new StringReader(props.toString()));
    InterpolationWarpScriptExtension ext = new InterpolationWarpScriptExtension();
    WarpScriptLib.register(ext);

    MemoryWarpScriptStack stack = new MemoryWarpScriptStack(null, null);
    stack.maxLimits();
    double[] xval = new double[] {1, 2};
    double[] yval = new double[] {1, 4};
    double[] zval = new double[] {1, 6};

    List<Double> xlist = new ArrayList<>();
    for (int i = 0; i < xval.length; i++) {
      xlist.add(xval[i]);
    }
    List<Double> ylist = new ArrayList<>();
    for (int i = 0; i < yval.length; i++) {
      ylist.add(yval[i]);
    }
    List<Double> zlist = new ArrayList<>();
    for (int i = 0; i < zval.length; i++) {
      zlist.add(zval[i]);
    }

    // Function values
    TrivariateFunction f = new TrivariateFunction() {
      @Override
      public double value(double x, double y, double z) {
        return 2 * x - 3 * y - 4 * z + 5;
      }
    };

    double[][][] fval = new double[xval.length][yval.length][zval.length];
    List<List<List<Double>>> fvalTensor = new ArrayList<List<List<Double>>>();

    for (int i = 0; i < xval.length; i++) {
      fvalTensor.add(new ArrayList<>());

      for (int j = 0; j < yval.length; j++) {
        fvalTensor.get(i).add(new ArrayList<Double>());

        for (int k = 0; k < zval.length; k++) {
          // tensor
          List<Double> l = fvalTensor.get(i).get(j);
          l.add(f.value(xval[i], yval[j], zval[k]));

          // 3d arrays
          fval[i][j][k] = f.value(xval[i], yval[j], zval[k]);
        }
      }
    }

    stack.push(xlist);
    stack.push(ylist);
    stack.push(zlist);
    stack.push(fvalTensor);

    ((WarpScriptStackFunction) ext.getFunctions().get(TRICUBICFIT.class.getSimpleName())).apply(stack);
    WarpScriptStackFunction function = (WarpScriptStackFunction) stack.pop();

    List<Number> input = new ArrayList<>();
    input.add(1.6);
    input.add(1);
    input.add(1);
    stack.push(input);
    function.apply(stack);
    stack.push(f.value(input.get(0).doubleValue(), input.get(1).doubleValue(), input.get(2).doubleValue()));

    System.out.println(stack.dump(100));
    //stack.execMulti("== ASSERT");

    //System.out.println(stack.dump(100));
  }

  /**
   * Test case for the microsphere interpolator
   */
  @Test
  public void microsphere_test() throws Exception {
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
    WarpConfig.setProperty(MICROSPHEREFIT.CONFIG_OR_CAPNAME_MAX_ELEMENTS, "100");
    stack.execMulti("{ 'elements' 100 'exponent' 1.1 } MICROSPHEREFIT");
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

  @Test
  public void startTestingWarp10Platform() throws Exception {

    final String HOME = "/home/jenx/Softwares/warp10-standalone/";
    //final String HOME = "/home/jenx/Softwares/warp10-dev/";

    //
    // Base configuration
    //

    String default_conf_folder = HOME + "etc/conf.d";
    List<String> conf = Files.walk(Paths.get(default_conf_folder)).map(x -> x.toString()).filter(f -> f.endsWith(".conf")).collect(Collectors.toList());

    //
    // Additional or overwriting configurations
    //

    String extraConfPath = HOME + "etc/conf.d/99-extra.conf";
    FileWriter fr = new FileWriter(new File(extraConfPath));
    fr.write("warp.timeunits = us\n");
    fr.close();
    conf.add(extraConfPath);

    //
    // Logging
    //

    System.setProperty("log4j.configuration", new File(HOME + "etc/log4j.properties").toURI().toString());

    //
    // Extensions and plugins from jars
    //

    List<File> jars = new ArrayList<File>();
    for (File f: new File(HOME + "lib").listFiles()) {
      if (f.getName().substring(f.getName().length() - 3).equals("jar")) {
        jars.add(f);
      }
    }

    URLClassLoader cl = (URLClassLoader) WarpScriptLib.class.getClassLoader();
    Method m = URLClassLoader.class.getDeclaredMethod("addURL", URL.class);
    m.setAccessible(true);

    for (int i = 0; i < jars.size(); i++) {
      URL url = jars.get(i).toURL();
      m.invoke(cl, url);
      System.out.println("Loading " + url.toString());
    }

    //
    // Start Warp 10
    //

    Thread t = new Thread(new Runnable() {
      @Override
      public void run() {
        try {
          Warp.main(conf.toArray(new String[conf.size()]));
        } catch (Exception e) {
          throw new RuntimeException(e);
        }
      }
    });
    t.start();

    //
    // Built-in extensions and plugins
    //

    WarpScriptLib.register(new InventoryWarpScriptExtension());
    //WarpScriptLib.register(new HttpWarpScriptExtension());

    while(null == Warp.getKeyStore()) {}
    WarpScriptLib.register(new TokenWarpScriptExtension(Warp.getKeyStore())); //null exception, must be loaded with config file
    System.out.println("Loaded Token extension");

    //
    // Extensions and plugins from the classpath
    //

    WarpScriptLib.register(new InterpolationWarpScriptExtension());

    //
    // Let Warp10 run
    //

    t.join();
  }

}