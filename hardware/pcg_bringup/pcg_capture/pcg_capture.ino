/*
 * pcg_capture.ino
 * MAX4466 fetal phonocardiography front-end: raw capture
 * ESP32 devkit, GPIO39 (ADC1_CH3, silkscreened "SN" = SENSOR_VN)
 * Streams CSV over serial: t_us,adc  @ 500 Hz
 *
 * DERIVED FROM toco_capture.ino WITH THREE CHANGES ONLY.
 * That sketch was measured at 500.000 Hz achieved with 9.4 us jitter std over
 * a 90 s characterisation run, and survived a 20-minute capture. Rewriting it
 * would discard that evidence, so the structure is left alone:
 *   - pin 35 -> 39, and the header string updated to match
 *   - settling delay 200 ms -> 2000 ms
 *   - comments
 * The timing loop, the signed-cast wraparound guard, the ADC configuration and
 * the line format are untouched.
 *
 * WHY 500 Hz
 * Not a compromise. The frozen pipeline's bandpass evaluates its upper corner
 * as min(200, 0.95*Nyquist); at 500 Hz that is 25-200 Hz, the IDENTICAL
 * passband used on 1000 Hz simfpcgdb data, and wider than the 25-158 Hz the
 * pipeline used on 333 Hz fpcgdb. 1 kHz would also not fit the serial link:
 * 15 bytes/line x 1000 = 15000 B/s against 11520 B/s available at 115200 8N1.
 * At 500 Hz it is 7500 B/s, ~65% utilisation, the figure the toco phase
 * measured as stable.
 *
 * WHY micros() AND NOT millis()
 * The AD8232 firmware timestamped with millis() at 1 ms resolution while
 * scheduling on micros(). Sub-millisecond jitter was invisible by
 * construction, and a "zero jitter" claim had to be withdrawn. PCG runs at
 * 4x the ECG rate, so the sample period here is 2000 us and millis() would
 * quantise it to two distinct values.
 *
 * NO SENTINEL COLUMN
 * The AD8232 firmware wrote -1 on lead-off. That value is FINITE, so it passed
 * straight through clean_nonfinite() as a huge negative spike instead of being
 * interpolated. A microphone has no lead-off condition to encode, so the
 * failure mode cannot recur here. If any fault flag is ever needed, it goes in
 * a SEPARATE COLUMN, never in the ADC value.
 *
 * TIMESTAMPS ARE RELATIVE TO t0
 * Keeps the field at 8 digits for captures under ~16 minutes, which holds the
 * line at 15 bytes. Ten-digit timestamps would push utilisation past 74%.
 * micros() wraps at ~71.6 minutes; the unsigned subtraction (t - t0) stays
 * correct across one wrap, but a capture longer than that would need handling.
 * All planned captures are 60-90 s.
 */

const int      MIC_PIN   = 39;                      // SN on the silkscreen.
                                                    // Single-channel erratum
                                                    // test showed no glitch
                                                    // signature vs G34; the
                                                    // MULTI-channel case is
                                                    // still untested.
const uint32_t SAMPLE_HZ = 500;
const uint32_t PERIOD_US = 1000000UL / SAMPLE_HZ;   // 2000 us

uint32_t t0;
uint32_t next_us;

void setup() {
  Serial.begin(115200);
  analogReadResolution(12);
  analogSetPinAttenuation(MIC_PIN, ADC_11db);

  // 2000 ms, up from toco's 200 ms. The MAX4466 output is AC-coupled around a
  // VCC/2 bias that has to charge on power-up. That time constant has NOT been
  // measured and no claim is made about it. The smoke runs showed elevated p2p
  // in the first 2-3 seconds, but the same pattern recurred mid-run, so it
  // could not be attributed to settling rather than room noise. The delay is
  // precautionary and free; the capture script trims a pre-registered lead-in
  // regardless.
  delay(2000);

  Serial.println("# pcg_capture max4466 pin=39 fs=500 att=11db");
  Serial.println("t_us,adc");
  t0      = micros();
  next_us = t0;
}

void loop() {
  next_us += PERIOD_US;
  while ((int32_t)(micros() - next_us) < 0) { }   // hold until the slot opens.
                                                  // Signed cast makes the
                                                  // comparison wraparound-safe.

  uint32_t t = micros();      // read AFTER the wait, so the timestamp records
                              // when the sample was actually taken, not when
                              // it was scheduled. Jitter stays visible.
  int v = analogRead(MIC_PIN);

  Serial.print(t - t0);
  Serial.print(',');
  Serial.println(v);
}
