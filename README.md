# North
Figure out which way's up from x-ray test shots. 

### Dependencies
 - [DIALS](https://dials.github.io)
 - [numpy](https://numpy.org/)
 - [pandas](https://pandas.pydata.org/)
 - [gspread](https://gspread.readthedocs.io/en/latest/)
 - [gspread-dataframe](https://pypi.org/project/gspread-dataframe/)
 - [oauth2client](https://oauth2client.readthedocs.io/en/latest/)

### Installation
  - `conda create -n north python=3.7`
  - `pip install gspread-dataframe`
  - `pip install --upgrade oauth2client`

### Contents
- scripts directory
	Contains useful dials shell scripts
- server.py
	Loops over a directory hierarchy and looks for new data to analyze and summarize
