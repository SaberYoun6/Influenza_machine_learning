#import TheBigKahuna2
from Tkinter import *
import subprocess as sub
class App:
    def __init__(self,master):

        frame = Frame(master)

        frame.pack()


        self.button = Button(
            frame, text='QUIT', fg='red', command=frame.quit
            )
        self.button.pack(side=RIGHT)
        
        self.run_the_big= Button(
            frame, text='Run', command=self.run_the_big
          )
        self.run_the_big.pack(side=LEFT)
        
	self.obtain_file= Button(
            frame, text="open file", command=self.obtain_file(f))## f equal files 
        self.obtain_file.pack(side=LEFT) 
       
        self.creating_file= Button(
            frame, text= "save file", command=self.creating_file(f))## f eqauls files
        self.creating_file.pack(side=LEFT)	
     
    def obtain_file(self,filen):
        if filen:
            print(filen.name)
    def creating_file(self,filename):
        if filename:
            print(filename.name)
    def run_the_big(self):
        print sub.call(["python","TheBigKahuna2.py","H1N11718.fasta","output.csv"])

root =Tk()

app=App(root)



root.mainloop() 
root.destroy()
