#import TheBigKahuna2
from Tkinter import *
from tkFileDialog import * 
import subprocess as sub
import csv
import os

class App:
    def __init__(self,master):

        frame = Frame(master)

        frame.pack()

        self.button = Button(
            frame, text='QUIT', fg='red', command=frame.quit
            )
        self.button.pack(side=RIGHT)
        
        self.run_the_big= Button(
            root, text='Run', command=self.run_the_big
          )
        self.run_the_big.pack(side=LEFT)
        
	self.obtain_file= Button(
            root, text="open file", command=self.obtain_file
          )## f equal files 
        self.obtain_file.pack(side=LEFT) 
       
        self.creating_file= Button(
            root, text= "save file", command=self.creating_file)## f eqauls files
        self.creating_file.pack(side=LEFT)


    def obtain_file(self):
        fToOpen = askopenfilename(filetypes=[("Fasta file","*.fasta"),("ALL files","*.*")])
        if fToOpen:
            try: 
                print("""here it comes: self.settings["template".set(fname)""")
            except:("Open source File", "Failed to read file\n'%s'" %fname)

        fileToOpen=open(fToOpen, 'r')
        return fileToOpen
        Label(root, fileToOpen.read()).pack()
        fToOpen.close()
    def creating_file(self):
        fToSave= asksaveasfile(mode='w',filetypes=[("comman separate values",'*.csv'),("All files","*.*")])
        if fToSave is None:
            return 
        with open(fToSave,'wb') as fp:
            fp.write('\ufeff')
            fileToSave=csv.writer(open(fp ,quoting=csv.QUOTE_MINIMAL))
            fileToSave.writerow(['virus type A strains', ",gc content","CpG ratio","GRAVY",["Amino"]*100])
        buffer =open(fToSave, 'r').read()
        print (buffer)
        f=open(fToSave, "w")
        f.write(buffer)
        f.close()

    def run_the_big(self):
        print sub.call(["python","TheBigKahuna2.py" ,"H3N21617.fasta", "output3.csv"])
            #self.obtain_file() , self.create_file()])

cwd= os.getcwd()
root =Tk()
root.geometry("500x500")

app=App(root)

root.mainloop() 
root.destroy()
