/*
 * pcg_smoke.ino
 * MAX4466 first-power diagnostic. ESP32 devkit, GPIO34 (ADC1_CH6).
 *
 * THIS IS NOT A CAPTURE. Nothing is timestamped, nothing is recorded, and no
 * claim in the notes may be sourced from it. Its only jobs are:
 *
 *   1. Confirm the board is alive and the DC bias sits where it should.
 *   2. Give a measured basis for setting the gain trimpot, which starts at an
 *      unknown factory position somewhere in the 25x-125x range. Setting gain
 *      by ear or by eye is guesswork; these statistics make it a measurement.
 *   3. Establish a quiet-room noise floor to compare later captures against.
 *
 * Prints one summary line per second at 115200 baud:
 *     sec,n,min,max,mean,std,p2p
 *
 * Sampling is paced from micros() at 500 Hz, matching the intended capture
 * rate, so the noise floor measured here is the noise floor the capture will
 * see. Statistics are accumulated on-device and only the summary is sent, so
 * the serial link is nowhere near its limit and cannot distort the timing.
 */

const int      MIC_PIN   = 39;                      // G34 on the silkscreen
                                                    // erratum, so anything odd
                                                    // here belongs to the mic,
                                                    // the gain, or the wiring.
const uint32_t SAMPLE_HZ = 500;
const uint32_t PERIOD_US = 1000000UL / SAMPLE_HZ;   // 2000 us

void setup() {
  Serial.begin(115200);
  analogReadResolution(12);
  analogSetPinAttenuation(MIC_PIN, ADC_11db);       // same as toco_capture
  delay(200);
  Serial.println("# pcg_smoke max4466 pin=39 fs=500 att=11db");
  Serial.println("# sec,n,min,max,mean,std,p2p");
}

void loop() {
  static uint32_t sec = 0;

  int      vmin = 4095, vmax = 0;
  double   sum = 0.0, sumsq = 0.0;
  uint32_t n = 0;

  uint32_t next_us = micros();
  for (uint32_t i = 0; i < SAMPLE_HZ; i++) {
    next_us += PERIOD_US;
    while ((int32_t)(micros() - next_us) < 0) { }   // signed cast: wraparound-safe
    int v = analogRead(MIC_PIN);
    if (v < vmin) vmin = v;
    if (v > vmax) vmax = v;
    sum   += v;
    sumsq += (double)v * (double)v;
    n++;
  }

  double mean = sum / n;
  double var  = sumsq / n - mean * mean;
  if (var < 0) var = 0;                             // guard rounding

  Serial.print(sec++);      Serial.print(',');
  Serial.print(n);          Serial.print(',');
  Serial.print(vmin);       Serial.print(',');
  Serial.print(vmax);       Serial.print(',');
  Serial.print(mean, 1);    Serial.print(',');
  Serial.print(sqrt(var), 2); Serial.print(',');
  Serial.println(vmax - vmin);
}
