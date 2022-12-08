use std::thread;
use std::sync::mpsc::channel;
use ndarray::{Array3, Array2};
use ndarray_npy::write_npy;

use crate::observer::{Observer, Photon, PhotonData, PHOTON_BATCH_SIZE, WIDTH, HEIGHT};
use crate::metrics::Metric;
use crate::source::AccretionDisk;

const PARALLELISM: bool = true;

enum ToThread {
    Photons([Option<Photon>; PHOTON_BATCH_SIZE]),
    Terminate(),
}

struct FromThread {
    results: Vec<PhotonData>,
    thread_index: usize,
}

pub struct Engine<O: Observer> {
    observer: O,
    optical_name: String,
    xray_name: String,
    optical: Array3<f64>,
    xray: Array3<f64>,
    counts: Array2<f64>,
    photon_count: usize,
}

impl<O: Observer> Engine<O> {
    pub fn new(observer: O, file_name: String) -> Self {
        Self {
            observer,
            optical_name: format!("../data/{}-optical.npy", file_name),
            xray_name: format!("../data/{}-xray.npy", file_name),
            optical: Array3::zeros((3, HEIGHT, WIDTH)),
            xray: Array3::zeros((3, HEIGHT, WIDTH)),
            counts: Array2::zeros((HEIGHT, WIDTH)),
            photon_count: 0,
        }
    }

    pub fn save(&self) {
        write_npy(&self.optical_name, &(&self.optical / &self.counts)).unwrap();
        write_npy(&self.xray_name, &(&self.xray / &self.counts)).unwrap();
    }

    pub fn run<M: Metric>(&mut self, max_iterations: usize, dtau: f64, metric: M, source: AccretionDisk) -> usize {
        let num_threads = if PARALLELISM {
            thread::available_parallelism().unwrap().get() - 1
        } else {
            1
        };
        println!("Using {} threads", num_threads);
        let mut threads = Vec::with_capacity(num_threads);
        let mut senders = Vec::with_capacity(num_threads);
        let (from_threads_sender, from_threads_receiver) = channel();
        let rng = fastrand::Rng::new();

        for thread_index in 0..num_threads {
            let (to_thread_sender, to_thread_receiver) = channel();
            let this_sender = from_threads_sender.clone();
            let this_source = source.clone();
            let this_rng = rng.clone();

            // Spawn the threads
            let t = thread::spawn(move || {
                while let Ok(msg) = to_thread_receiver.recv() {
                    // Handle message
                    let mut results = Vec::with_capacity(PHOTON_BATCH_SIZE);
                    match msg {
                        ToThread::Photons(photons) => {
                            for photon in photons {
                                results.push(match photon {
                                    Some(p) => p.run(max_iterations, dtau, &metric, &this_source, &this_rng),
                                    None => break,
                                });
                            }
                        },
                        ToThread::Terminate() => {
                            return;
                        }
                    };
                    
                    // Commit results
                    this_sender.send(FromThread {
                        thread_index, 
                        results,
                    }).unwrap();
                }
            });

            to_thread_sender.send(ToThread::Photons(self.observer.next_photons())).unwrap();

            threads.push(t);
            senders.push(to_thread_sender);
        }

        // Main loop
        let mut last_percentage = 0;
        while let Ok(msg) = from_threads_receiver.recv() {
            // Update the observer with new info
            let percentage = self.observer.update(&msg.results);
            if percentage > last_percentage {
                if percentage % 10 == 0 {
                    println!("{}%", percentage);
                }
                last_percentage = percentage;
            }
            if percentage >= 100 {
                break;
            }
            senders[msg.thread_index].send(ToThread::Photons(self.observer.next_photons())).unwrap(); // Start new results

            self.add(msg.results);
        }
        
        // Clean up
        for sender in senders {
            sender.send(ToThread::Terminate()).unwrap();
        }
        for thread in threads {
            thread.join().unwrap();
        }
        // Consume the very last photons
        while let Ok(msg) = from_threads_receiver.try_recv() {
            // Update the observer with new info
            self.observer.update(&msg.results);
            self.add(msg.results);
        }

        self.save();

        self.photon_count
    }

    fn add(&mut self, results: Vec<PhotonData>) {
        self.photon_count += results.len();
        for p in results {
            self.optical[(0, p.pixel.0,p.pixel.1)] += p.optical_color.0;
            self.optical[(1, p.pixel.0,p.pixel.1)] += p.optical_color.1;
            self.optical[(2, p.pixel.0,p.pixel.1)] += p.optical_color.2;
            self.xray[(0, p.pixel.0,p.pixel.1)] += p.xray_color.0;
            self.xray[(1, p.pixel.0,p.pixel.1)] += p.xray_color.1;
            self.xray[(2, p.pixel.0,p.pixel.1)] += p.xray_color.2;
            self.counts[(p.pixel.0,p.pixel.1)] += 1.0;
        }
    }
}