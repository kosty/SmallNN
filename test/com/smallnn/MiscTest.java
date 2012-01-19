/**
 *  This file is part of SmallNN, a small neural network implementation
 *  Copyright (C) 2011, 2012 Arsen Kostenko <arsen.kostenko@gmail.com>
 *     
 *  SmallNN is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SmallNN is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with SmallNN.  If not, see <http://www.gnu.org/licenses/>.
 */
package com.smallnn;

import static java.util.Arrays.binarySearch;
import static junit.framework.Assert.assertEquals;

import java.util.Arrays;

import org.junit.Test;

public class MiscTest {

    @Test
    public void test() {
        double[] rates = {1.3, 0.9, 0.6, 0.3, 0.09, 0.06, 0.03, 0.009, 0.006, 0.003, 0.0009, 0.0006, 0.0003};
        Arrays.sort(rates);
        int idx = binarySearch(rates, 0.6d);
        assertEquals(2, idx);
    }

}
