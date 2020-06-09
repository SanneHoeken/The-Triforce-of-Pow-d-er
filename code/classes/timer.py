import time

class Timer():

    def __init__(self):
        self.start_time = None
        self.passed_time = None

    def start(self):
        self.start_time = time.perf_counter()

    def stop(self):
        self.passed_time = round(time.perf_counter() - self.start_time, 4)
        self.start_time = None

    def get_time(self):
        return self.passed_time