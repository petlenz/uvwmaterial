.. Copyright (c) 2022, Peter Lenz

   Distributed under the terms of the BSD 3-Clause License.

   The full license is in the file LICENSE, distributed with this software.
   
Getting started
===============

This guide explains how to get started with ``uvwmaterial`` after installing it 
using the methods described in the installation section.



First example
-------------
.. code::

    #include <iostream>
    #include <tmech/tmech.h>
    #include <uvwmaterial/uvwmaterial.h>
    
    int main()
    {

        return 0;
    }

This example ...




Compiling the first example
---------------------------

``uvwmaterial`` is a header-only library. You have to tell the compiler where to find ``uvwmaterial``.
For example with g++, use the ``-I`` option to achieve this. Assuming the first example code is 
located in ``first_example.cpp``, the compilation command is:

.. code:: bash

    g++ -std=c++17 -I /path/to/tmech/ -I /path/to/uvwmaterial/ first_example.cpp -o first_example


Starting the program, produces the following output:

.. code::


Building with cmake
-------------------

A better way for building programs using ``uvwmaterial`` is to use `cmake`. 
Assuming the following folder structure:

.. code:: bash

    first_example
       |- src
       |   |- first_example.cpp
       |- CMakeLists.txt

The following minimal ``CMakeLists.txt`` is enough to build the first example:

.. code:: cmake

    cmake_minimum_required(VERSION 3.1)
    set(CMAKE_CXX_STANDARD 17)
    project(first_example)
    
    find_package(tmech REQUIRED)
    find_package(uvwmaterial REQUIRED)

    add_executable(first_example src/first_example.cpp)

    if(MSVC)
        set(CMAKE_EXE_LINKER_FLAGS /MANIFEST:NO)
    endif()

    target_link_libraries(first_example tmech uvwmaterial)


`cmake` has to know where to find the headers, this is done through the ``CMAKE_INSTALL_PREFIX``
variable. Note that ``CMAKE_INSTALL_PREFIX`` is usually the path to a folder containing the following
subfolders: ``include``, ``lib`` and ``bin``, so you don't have to pass any additional option for linking.
Examples of valid values for ``CMAKE_INSTALL_PREFIX`` on Unix platforms are ``/usr/local``, ``/opt``.

The following commands create a directory for building, builds
the first example with cmake and then runs the program:

.. code:: bash

    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=your_prefix ..
    make
    ./first_program


Second example
--------------

.. code::

    #include <iostream>
    #include <tmech/tmech.h>

    int main()
    {
        return 0;
    }

When compiled and run, this produces the following output:

.. code::


Third example
-------------

.. code::

    #include <iostream>
    #include <tmech/tmech.h>

    int main()
    {
        return 0;
    }

Outputs:

.. code::

