The program generates a Mathieu Stability Diagram for a specific mass and mass tolerance.
Matplotlib is used for plot generation and PyQT for the window.


TODO:

General:
    - More commenting

File system:
    - Read exp setup from file?
    - Read/save latest specific setup from/to file

Specific functionality:
    - Smoother tuning / non linear DC tuning?
    - Check stability of masses - with color coordination

Code cleanup
    - remove functionality overlap in methods
    - Maybe? Port matplotlib code to PyQT plotting
    - Split code into frontend plotting/UI and backend mathieu calculations?  