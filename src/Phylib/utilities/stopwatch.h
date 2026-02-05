//
//  StopWatch.h
//
// Structure to help with timing
//  Created by David Bryant on 16/07/2025.
//
// Changed clock type 5/8/25
//

#include <chrono>

/**
 Stopwatch class
 
 start()    Set the start time to now.
 stop()    Figure out the time since start was pressed and add this to the elapsed time
 reset()   Set the elapsed time to zero
get()  If running, add the current time to elapsed time and output that, otherwise just output the elapsed time. (in milliseconds)
 */
class Stopwatch {
public:
    Stopwatch() : running(false), elapsed(0) {}

    void start() {
        startTime = std::chrono::steady_clock::now();
        running = true;
    }

    void stop() {
        if (running) {
            auto endTime = std::chrono::steady_clock::now();
            elapsed += (double) std::chrono::duration_cast<std::chrono::microseconds>(endTime-startTime).count()/1000.0;
            running = false;
        }
    }

    void reset() {
        elapsed = 0;
        running = false;
    }

    double get() const {
        if (running) {
            auto now = std::chrono::steady_clock::now();
            return elapsed + (double)std::chrono::duration_cast<std::chrono::microseconds>(now-startTime).count()/1000.0;
        } else {
            return elapsed;
        }
    }

private:
    bool running;
    double elapsed; // milliseconds
    std::chrono::steady_clock::time_point startTime;
};

