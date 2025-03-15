.include "m8515def.inc"

.def reg_led = r20
.def temp    = r19
.def direction = r22
.def job = r21          

.equ START = 0
.equ STOP = 1

.org $000
    rjmp INIT             ; RESET
    reti                  ; INT0 
    rjmp START_ISR        ; INT1 - запуск
    rjmp STOP_ISR         ; INT2 - остановка

INIT:
    ldi temp, high(RAMEND)
    out SPH, temp
    ldi temp, low(RAMEND)
    out SPL, temp
=
   
    ser temp
    out DDRB, temp     
    clr temp
    out DDRD, temp      
    out DDRE, temp       

    ldi temp, (1<<PD3)
    out PORTD, temp
    ldi temp, (1<<PE0)
    out PORTE, temp

    
    ldi temp, (1<<INT1)|(1<<INT2)
    out GICR, temp

    clr temp
    out MCUCR, temp       ; INT1 
    out EMCUCR, temp      ; INT2 

    sei                

    ldi reg_led, 0b00000011
    ldi direction, 0    
    clr job              

LOOP:
    sbrs job, 0       
    rjmp LOOP            

    com reg_led
    out PORTB, reg_led
    com reg_led         

    rcall DELAY_MS       

    tst direction
    breq MOVE_RIGHT
    rjmp MOVE_LEFT

MOVE_RIGHT:
    lsl reg_led
    cpi reg_led, 0b11000000
    brne LOOP
    ldi direction, 1
    rjmp LOOP

MOVE_LEFT:
    lsr reg_led
    cpi reg_led, 0b00000011
    brne LOOP
    clr direction
    rjmp LOOP

START_ISR:
    ldi job, 1
wait_start_release:
    sbic PIND, PD2
    rjmp wait_start_release
    reti

; INT1 - остановка
STOP_ISR:
    clr job
wait_stop_release:
    sbic PIND, PD3
    rjmp wait_stop_release
    reti


DELAY_MS:
    ldi r16,9
del3:
    ldi r17,202
del2:
    ldi r18,221
del1:
    dec r18
    brne del1
    dec r17
    brne del2
    dec r16
    brne del3
    ret
