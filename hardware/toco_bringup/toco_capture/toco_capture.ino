/*
 * toco_capture.ino
 * FSR402 uterine activity front-end: raw characterisation capture
 * ESP32 devkit, GPIO35 (ADC1_CH7), 10k pull-down divider
 * Streams CSV over serial: t_us,adc  @ 500 Hz
 */

const int      FSR_PIN   = 35;
const uint32_t SAMPLE_HZ = 500;
const uint32_t PERIOD_US = 1000000UL / SAMPLE_HZ;   // 2000 us

uint32_t t0;
uint32_t next_us;

void setup() {
  Serial.begin(115200);
  analogReadResolution(12);
  analogSetPinAttenuation(FSR_PIN, ADC_11db);
  delay(200);
  Serial.println("# toco_capture fsr402 pin=35 fs=500 rfixed=10k");
  Serial.println("t_us,adc");
  t0      = micros();
  next_us = t0;
}

void loop() {
  next_us += PERIOD_US;
  while ((int32_t)(micros() - next_us) < 0) { }   // hold until the slot opens

  uint32_t t = micros();
  int v = analogRead(FSR_PIN);

  Serial.print(t - t0);
  Serial.print(',');
  Serial.println(v);
}