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
import io.warp10.script.ext.inventory.InventoryWarpScriptExtension;
import io.warp10.script.ext.token.TokenWarpScriptExtension;
import io.warp10.standalone.Warp;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.junit.Assert;
import org.junit.Ignore;
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
import java.util.List;
import java.util.stream.Collectors;

public class InterpolationWarpScriptExtensionTest {

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