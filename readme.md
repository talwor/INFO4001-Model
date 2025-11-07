Python 3.10 is used for this project, other versions should be fine but 3.10 were used to develop this project

PACKAGES REQUIRED TO BE INSTALLED VIA PIP

pip install networkx
pip install matplotlib
pip install random


If on mac, you might encounter module not found error like I have, 
I have navigated this problem by using a virtual environment and installing the dependencies in said network.
I have provided a copy and paste below for a fast set up
Run each line independently one after each other

""""
1. /Library/Frameworks/Python.framework/Versions/3.10/bin/python3 -m venv venv 
2. source venv/bin/activate
3. pip install networkx
""" 


There are some control_tests also made, to run ONLY a singular disease or even both, with configurable days and seeds
Below is what I used for my thesis, but the numbers are interchangeable.

python control_tests.py --mode hiv_only --days 730 --seeds 0-20
python control_tests.py --mode flu_only --days 730 --seeds 0-20
python control_tests.py --mode both --days 730 --batch --seeds 0-20