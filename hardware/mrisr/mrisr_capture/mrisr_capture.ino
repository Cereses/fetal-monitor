/*
 * mrisr_capture.ino
 * Multi-rate acquisition: MAX4466 (PCG) + FSR402 (UC) on ONE time base.
 * ESP32 devkit. GPIO39 = PCG (silkscreen "SN"), GPIO35 = UC.
 *
 * WHY A SHARED TIME BASE IS THE POINT
 * A fetal heart rate deceleration is interpreted by its PHASE relative to the
 * contraction. Early decelerations (nadir at the contraction peak) are benign
 * head compression; late decelerations (nadir AFTER the peak) indicate
 * uteroplacental insufficiency. Same shape, opposite clinical meaning, and
 * only the timing separates them. Two independent clocks make that phase
 * unrecoverable, so FHR and contraction timing MUST share a time base or the
 * pairing -- which is what cardiotocography IS -- means nothing.
 *
 * DESIGN
 * One 500 Hz schedule, derived exactly as in pcg_capture.ino (measured at
 * 500.0000 Hz with 0.6 us jitter and zero gaps over 225k samples). Every tick
 * reads BOTH channels. PCG is emitted immediately. UC accumulates and is
 * emitted every 125th tick as a sum.
 *
 * WHY 125-SAMPLE AVERAGING AND NOT EVERY-125th SELECTION
 * This is a compatibility requirement, not a preference. 04_detect.py records
 * that device UC data is "float-valued (means of 125 raw samples), so
 * exact-equality flat runs essentially never occur" -- which makes
 * uc_detector's flat_run_mask INERT on device data, a property the toco notes
 * documented honestly rather than hiding.
 *
 * Taking every 125th sample instead would restore INTEGER values, silently
 * REACTIVATING that mask. Same frozen detector, different behaviour, no error
 * message. Averaging reproduces exactly what the frozen detector was validated
 * against. It is also a real (if weak) anti-alias filter: a 125-tap boxcar has
 * nulls at multiples of 4 Hz.
 *
 * THE SUM IS SENT, NOT THE MEAN
 * The host divides by 125.0. Avoids on-device float formatting and its
 * rounding, and reproduces the toco phase's float means bit-for-bit.
 * Max value 125 * 4095 = 511,875, comfortably inside uint32_t.
 *
 * UC TIMESTAMP IS THE WINDOW CENTRE
 * (t_first + t_last) / 2. A mean represents the centre of its window. Tagging
 * it with the last tick would introduce a systematic 125 ms lag between
 * channels -- and cross-channel phase is the entire reason this firmware
 * exists.
 *
 * SERIAL BUDGET
 *   PCG  500 lines/s x 15 bytes = 7,500 B/s
 *   UC     4 lines/s x 19 bytes =    76 B/s
 *                                 ---------
 *                                 7,576 B/s = 65.8% of 11,520 B/s
 * 0.1 points above pcg_capture.ino, which ran clean across 225k samples.
 *
 * LINE FORMATS
 *   PCG:  t_us,adc        (starts with a digit)
 *   UC:   U,t_us,sum      (starts with 'U')
 * No sentinel values anywhere. The AD8232's -1 lead-off marker was finite and
 * passed straight through clean_nonfinite() as a huge negative spike; any
 * future fault flag goes in its own field, never in a sample value.
 *
 * WHAT THIS FIRMWARE IS FOR MEASURING
 * Two reads per tick on alternating ADC1 channels is exactly the condition
 * notes.md 6 could NOT clear: the SENSOR_VP/VN erratum is most reported under
 * channel alternation, and that test used single-channel continuous sampling.
 * Cross-channel settling is the other unknown -- the ESP32 has one SAR
 * converter behind a multiplexer, and reading too soon after a channel switch
 * can return a value contaminated by the previous channel. Both are measured
 * against existing single-channel baselines, not assumed.
 */

const int      PCG_PIN   = 39;                      // "SN" on the silkscreen
const int      UC_PIN    = 35;
const uint32_t SAMPLE_HZ = 500;
const uint32_t PERIOD_US = 1000000UL / SAMPLE_HZ;   // 2000 us
const uint32_t UC_DECIM  = 125;                     // 500 / 125 = 4 Hz

uint32_t t0;
uint32_t next_us;

uint32_t uc_sum   = 0;
uint32_t uc_n     = 0;
uint32_t uc_first = 0;
uint32_t uc_last  = 0;

void setup() {
  Serial.begin(115200);
  analogReadResolution(12);
  analogSetPinAttenuation(PCG_PIN, ADC_11db);
  analogSetPinAttenuation(UC_PIN,  ADC_11db);

  delay(2000);          // MAX4466 output is AC-coupled around a VCC/2 bias
                        // that must charge. Precautionary; the capture script
                        // trims a pre-registered lead-in regardless.

  Serial.println("# mrisr_capture pcg=39@500 uc=35@4 decim=125 "
                 "fs=500 att=11db uc=mean_of_125_as_sum");
  Serial.println("# pcg: t_us,adc   uc: U,t_us_centre,sum");
  t0      = micros();
  next_us = t0;
}

void loop() {
  next_us += PERIOD_US;
  while ((int32_t)(micros() - next_us) < 0) { }   // signed cast: wraparound-safe

  uint32_t t = micros();          // read AFTER the wait, so the timestamp is
                                  // when the sample was TAKEN, not scheduled.
                                  // Jitter stays visible.

  int pcg = analogRead(PCG_PIN);
  int uc  = analogRead(UC_PIN);   // channel switch happens here -- the thing
                                  // being measured

  Serial.print(t - t0);
  Serial.print(',');
  Serial.println(pcg);

  if (uc_n == 0) uc_first = t;
  uc_last = t;
  uc_sum += (uint32_t)uc;
  uc_n++;

  if (uc_n >= UC_DECIM) {
    Serial.print("U,");
    Serial.print(uc_first + (uc_last - uc_first) / 2 - t0);   // window centre
    Serial.print(',');
    Serial.println(uc_sum);
    uc_sum = 0;
    uc_n   = 0;
  }
}
