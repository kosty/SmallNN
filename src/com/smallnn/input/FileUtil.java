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
package com.smallnn.input;

import static org.apache.commons.lang.SystemUtils.getUserHome;

import java.io.File;
import java.io.FileFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class FileUtil {

    public static List<File> find(File file) {
        return find(file, null);
    }
    
    public static final FileFilter pngFilter  = new FileFilter() {

        @Override
        public boolean accept(File f) {
            return f.isFile() && f.getName().toLowerCase().endsWith("png");
        }
        
    };
    
    public static final FileFilter dirFilter  = new FileFilter() {

        @Override
        public boolean accept(File f) {
            return f.isDirectory();
        }
        
    };
    
    public static List<File> find(File file, FileFilter filter) {
        List<File> result = new ArrayList<File>();
        LinkedList<File> stack = new LinkedList<File>();
        stack.push(file);
        while (!stack.isEmpty()) {
            File f = stack.pop();
            if (filter == null || filter.accept(f)) {
                result.add(f);
            }

            if (f.isDirectory() && f.exists()) {
                stack.addAll(Arrays.asList(f.listFiles()));
            }
        }
        return result;
    }

    public static File tildeExpand(String path) {
        if (path.startsWith("~")) {
            path = path.replaceFirst("~", getUserHome().getAbsolutePath());
        }
        return new File(path);
    }
}
