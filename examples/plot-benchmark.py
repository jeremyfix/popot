
import sys
import json
import glob
import matplotlib.pyplot as plt
import numpy as np
import PyQt4.QtCore as pyqtcore
import PyQt4.QtGui as pyqtgui


class MyTable(pyqtgui.QTableWidget):
    
    def __init__(self, json_files, qApp, *args):
        pyqtgui.QTableWidget.__init__(self, len(json_files), 5)
        self.qApp = qApp
        horHeaders = ['Use?', 'Problem', 'Dimension', 'Algorithm', 'File']
        self.setHorizontalHeaderLabels(horHeaders)

        for i, f in enumerate(json_files):
            fp = open(f)
            results = json.load(fp)
            algo = results['algorithm']
            problem = results['problem']

            chkBoxItem = pyqtgui.QTableWidgetItem()
            chkBoxItem.setFlags(pyqtcore.Qt.ItemIsUserCheckable | pyqtcore.Qt.ItemIsEnabled)
            chkBoxItem.setCheckState(pyqtcore.Qt.Unchecked)       
            self.setItem(i,0,chkBoxItem)


            row_elements = [problem['name'], str(problem['dimension']), algo, f]
            for j, l in enumerate(row_elements):
                newitem = pyqtgui.QTableWidgetItem(l)
                newitem.setFlags(pyqtcore.Qt.NoItemFlags)
                self.setItem(i, j+1, newitem)
            fp.close()


        self.resizeColumnsToContents()
        self.resizeRowsToContents()
        self.setSortingEnabled(1)
        self.selected_files = []

    def buttonClicked(self):
        # We go through the table and list all the
        # json files that have been selected
        self.selected_files = []
        Nrows = self.rowCount()
        for i in range(Nrows):
            item = self.item(i, 0)
            if(item.checkState() == pyqtcore.Qt.Checked):
                self.selected_files.append(self.item(i, 4).text())
        self.selected_files = map(str, self.selected_files)
        self.qApp.quit()


json_files = []

if(len(sys.argv) == 1):
    # We list all the .json files in the current directory
    json_files = glob.glob('./*.json')
    use_selector = True
else:
    json_files = sys.argv[1:]
    use_selector = False


# We propose a selector where we list all the json files in the current directory
# you select the ones to plot 
if(use_selector):
    print("Using an interface to select the files to plot")
    print(json_files)
    
    app = pyqtgui.QApplication(sys.argv)
    frame	 = pyqtgui.QFrame() 
    vBoxLayout	 = pyqtgui.QVBoxLayout()
    table        = MyTable(json_files, app)
    button	 = pyqtgui.QPushButton("Plot the results")
    button.clicked.connect(table.buttonClicked)
    vBoxLayout.addWidget(table)
    vBoxLayout.addWidget(button)	

    frame.resize(640,480)
    frame.setLayout(vBoxLayout)
    frame.setWindowTitle("JSON result files selection")
    frame.updateGeometry()
    frame.show()	
    app.exec_()
    frame.hide()

    # The window has been closed, we recover the data files to plot
    json_files = table.selected_files



# Now we can plot all the elements within json_files
for f in json_files:
    fp = open(f)
    results = json.load(fp)

    error = np.array(results['mean'])
    std = np.array(results['std'])
    FE = np.array(results['FE'])
    algo = results['algorithm']
    problem = results['problem']
    title = "%s, d = %i" % (problem['name'], problem['dimension'])

    line, = plt.plot(np.log(FE)/np.log(10.), error, label=algo)
    plt.plot(np.log(FE)/np.log(10.), error+std,  '--', color=line.get_color())
    plt.plot(np.log(FE)/np.log(10.), error-std,  '--',color=line.get_color())


plt.ylabel(r'Mean function value')
plt.xlabel(r'Function evaluations (log)')

plt.legend()


plt.show()
