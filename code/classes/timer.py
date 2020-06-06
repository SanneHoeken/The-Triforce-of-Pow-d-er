import time

class Timer():

    def __init__(self):
        self.start_time = None

    def start(self):
        self.start_time = time.perf_counter()

    def stop(self):
        passed_time = time.perf_counter() - self.start_time
        self.start_time = None
        print(f"Time: {passed_time: 0.4f} seconds")