# Reaction Diffusion with Sverchok SNLite node in Blender

## Description
This project includes some examples of Reaction Diffusion (RD) simulation with SNLite node in Blender Sverchok add-on. 

This project includes:
- Python script that simulates 2D RD with C code.
- Python script that simulates 3D RD with C code.
- Python script that simulates 3D RD "without" C code.
- Python script that simulates RD on mesh with C code.
- Python script that simulates RD on mesh "without" C code.

Python scripts "with C code" uses these C codes as C dynamic libraries. So you have to compile these C codes to generate C dynamic libraries in advance. Though you have to do so, these scripts runs very fast and achieves high performance.

Scripts "without C code" can run without C library. But the processing speed and performance are relatively low.

This project idea is originated from [@zeffi's gist](https://gist.github.com/zeffii/9e156f0d37977fd1b0ca3c65d0ddc611) and [nortikin/sverchok#1734](https://github.com/nortikin/sverchok/issues/1734) .

## Usage

- Compile C code and make C dynamic library.
    
    ```console
    :: Windows
    cl.exe /D_USRDLL /D_WINDLL Basic_RD_xd.c /MT /link /DLL /OUT:libBasic_RD_xd.dll
    ```
    ```sh
    # macOS
    gcc -dynamiclib -o ./libBasic_RD_xd.dylib ./Basic_RD_xd.c
    ```
    ```sh
    # Linux
    gcc -c -fPIC Basic_RD_xd.c -o Basic_RD_xd.o
    gcc Basic_RD_xd.o -shared -o libBasic_RD_xd.so
    ```

- Launch Blender and open Sverchok node editor.
- Add Script Node Lite (SNLite) node.
- Open Text editor in Blender.
- Copy and paste one of the Python script to the text buffer.
- Modify "load_library" arguments in the script to load your library (if you use "with C code" script).
- Put the the text buffer name to the node and click the right button.

## Requirements
* Blender 2.8 (or later)
* sverchok add-on 0.6 (or later)

## Author
asahidari

## Licence
[GPL 3](https://www.gnu.org/licenses/quick-guide-gplv3.html)
