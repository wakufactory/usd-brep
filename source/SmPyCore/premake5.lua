-- SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
-- SPDX-License-Identifier: Apache-2.0

-- Small USD-free Python surface for local solid modeling and mesh export.
project "_smlib_core"
    kind "SharedLib"
    location (workspace_dir.."/%{prj.name}")
    targetprefix ""
    files { "src/**.cpp" }
    includedirs {
        root.."/source/SM_API/inc",
        root.."/source/SMLib/inc",
    }
    externalincludedirs { root.."/_build/target-deps/pybind11/include" }
    links { "SMLib", "SM_API" }

    filter { "system:windows" }
        externalincludedirs { root.."/_build/target-deps/python/include" }
        libdirs { root.."/_build/target-deps/python/libs" }
        targetextension ".pyd"
    filter { "system:linux", "platforms:x86_64" }
        externalincludedirs { root.."/_build/target-deps/python/include/python3.12" }
        libdirs { root.."/_build/target-deps/python/lib" }
        links { "python3.12" }
        targetextension ".cpython-312-x86_64-linux-gnu.so"
        buildoptions { "-fPIC", "-fvisibility=hidden" }
    filter { "system:linux", "platforms:aarch64" }
        externalincludedirs { root.."/_build/target-deps/python/include/python3.12" }
        libdirs { root.."/_build/target-deps/python/lib" }
        links { "python3.12" }
        targetextension ".cpython-312-aarch64-linux-gnu.so"
        buildoptions { "-fPIC", "-fvisibility=hidden" }
    filter { "system:macosx" }
        use_python()
        targetextension ".cpython-312-darwin.so"
        buildoptions (concat_arrays({"-fPIC", "-fvisibility=hidden", "-Wno-deprecated-literal-operator"}, macos_common_buildoptions()))
        linkoptions (concat_arrays({"-rpath @loader_path"}, macos_common_linkoptions()))
    filter {}
