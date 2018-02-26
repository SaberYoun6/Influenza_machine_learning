import multiprocessing as mpi
import time
import sys
import 

def worker():
    name = mpi.current_process().name
    print name , "Starting"
    time.sleep(60)
    print name, "Exiting"
def my_service():
    name = mpi.current_process().name
    print name, 'starting'
    time.sleep(50)
    print name, 'exiting'

def main():
    service = mpi.Process(name='my_service', target= my_service)
    work_1= mpi.Process(name='worker 1',target=worker)
    work_2= mpi.Process(target=worker)
    pool=mpi.Pool(4)
    work_1.start()
    work_2.start()
    service.start()
main()
