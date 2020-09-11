import RPi.GPIO as GPIO
from time import sleep             # lets us have a delay
GPIO.setmode(GPIO.BCM)             # choose BCM or BOARD
GPIO.setup(16, GPIO.OUT)           # set GPIO16 as an output

fake_pps = 16

while True:
    GPIO.output(fake_pps, 1)
    print("GPIO16 went on")
    sleep(0.1)
    GPIO.output(fake_pps, 0)
    print("GPIO16 went off")
    sleep(0.9)

