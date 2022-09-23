.. Copyright (c) 2022, Peter Lenz

   Distributed under the terms of the BSD 3-Clause License.

   The full license is in the file LICENSE, distributed with this software.
   
Install
=======

Although ``uvwmaterial`` is a header-only library, we provide standardized methods to install it with cmake.

Besides the ``uvwmaterial`` headers, all these methods place the ``cmake`` project
configuration file in the right location so that third-party projects can use
cmake's ``find_package`` to locate tmech headers.


From source with cmake
----------------------

You can also install ``uvwmaterial``` from source with cmake. On Unix platforms, from the
source directory:

.. code::

    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=path_to_prefix ..
    make install

On Windows platforms, from the source directory:

.. code::

    mkdir build
    cd build
    cmake -G "NMake Makefiles" -DCMAKE_INSTALL_PREFIX=path_to_prefix ..
    nmake
    nmake install

``path_to_prefix`` is the absolute path to the folder where cmake searches for
dependencies and installs libraries. ``uvwmaterial`` installation from cmake assumes
this folder contains ``include`` and ``lib`` subfolders.


Including ``uvwmaterial`` in your project
---------------------------------

After installing the library with cmake, you can add ``uvwmaterial`` to your project using cmake:

.. code::

    find_package(uvwmaterial REQUIRED)
    target_include_directories(your_target PUBLIC ${uvwmaterial_INCLUDE_DIR})
    target_link_libraries(your_target PUBLIC uvwmaterial)
