//
//   Copyright 2022 - 2023  SenX S.A.S.
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

import java.util.HashMap;
import java.util.Map;

import io.warp10.WarpConfig;
import io.warp10.script.WarpScriptStack;
import io.warp10.warp.sdk.Capabilities;
import io.warp10.warp.sdk.WarpScriptExtension;

public class InterpolationWarpScriptExtension extends WarpScriptExtension {

  //
  // Functions
  //

  public static final String INTERPOLATOR_ND_MICROSPHERE = "INTERPOLATOR.ND.MICROSPHERE";
  public static final String INTERPOLATOR_ND_SMICROSPHERE = "INTERPOLATOR.ND.SMICROSPHERE"; // seeded
  public static final String INTERPOLATOR_2D_BICUBIC = "INTERPOLATOR.2D.BICUBIC";
  public static final String INTERPOLATOR_3D_TRICUBIC = "INTERPOLATOR.3D.TRICUBIC";
  public static final String INTERPOLATOR_1D_LINEAR = "INTERPOLATOR.1D.LINEAR";
  public static final String INTERPOLATOR_1D_SPLINE = "INTERPOLATOR.1D.SPLINE";
  public static final String INTERPOLATOR_1D_AKIMA = "INTERPOLATOR.1D.AKIMA";
  public static final String INTERPOLATOR_TICK_LINEAR = "mapper.interpolator.linear";
  public static final String INTERPOLATOR_TICK_SPLINE = "mapper.interpolator.spline";
  public static final String INTERPOLATOR_TICK_AKIMA = "mapper.interpolator.akima";

  private static final Map<String, Object> functions;

  static {
    functions = new HashMap<String, Object>();

    functions.put(INTERPOLATOR_ND_MICROSPHERE, new InterpolatorMicrosphere(INTERPOLATOR_ND_MICROSPHERE, false));
    functions.put(INTERPOLATOR_ND_SMICROSPHERE, new InterpolatorMicrosphere(INTERPOLATOR_ND_SMICROSPHERE, true));
    functions.put(INTERPOLATOR_2D_BICUBIC, new InterpolatorBicubic(INTERPOLATOR_2D_BICUBIC));
    functions.put(INTERPOLATOR_3D_TRICUBIC, new InterpolatorTricubic(INTERPOLATOR_3D_TRICUBIC));
    functions.put(INTERPOLATOR_1D_LINEAR, new InterpolatorUnivariate(INTERPOLATOR_1D_LINEAR, InterpolatorUnivariate.TYPE.LINEAR));
    functions.put(INTERPOLATOR_1D_SPLINE, new InterpolatorUnivariate(INTERPOLATOR_1D_SPLINE, InterpolatorUnivariate.TYPE.SPLINE));
    functions.put(INTERPOLATOR_1D_AKIMA, new InterpolatorUnivariate(INTERPOLATOR_1D_AKIMA, InterpolatorUnivariate.TYPE.AKIMA));
    functions.put(INTERPOLATOR_TICK_LINEAR, new InterpolatorOnTicks(INTERPOLATOR_TICK_LINEAR, InterpolatorOnTicks.TYPE.LINEAR));
    functions.put(INTERPOLATOR_TICK_SPLINE, new InterpolatorOnTicks(INTERPOLATOR_TICK_SPLINE, InterpolatorOnTicks.TYPE.SPLINE));
    functions.put(INTERPOLATOR_TICK_AKIMA, new InterpolatorOnTicks(INTERPOLATOR_TICK_AKIMA, InterpolatorOnTicks.TYPE.AKIMA));
  }

  //
  // Function special configuration
  //

  static int getIntLimitValue(WarpScriptStack stack, String limitName, int defaultValue) {

    // this is the default limit value
    int val = defaultValue;

    // if the config exists with the limit name, it is used to determine the value of the limit
    if (null != WarpConfig.getProperty(limitName)) {
      val = Integer.parseInt(WarpConfig.getProperty(limitName));
    }

    // if the capability with the limit name is present within the stack, its value supersedes the default or configured limit
    String capValue = Capabilities.get(stack, limitName);
    if (null != capValue) {
      val = Integer.parseInt(capValue);
    }

    return val;
  }

  @Override
  public Map<String, Object> getFunctions() {
    return functions;
  }
}
