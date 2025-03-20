.include "m8515def.inc"
.def temp = r16
.equ led = 0     
.equ sw0 = 0         ; ������ SW0 � ������������� �� PA0
.equ sw1 = 1         ; ������ SW1 � ������������� �� PA1


.org 0x000
    rjmp INIT           
   	nop     
    nop 
    nop 
   	nop     
    nop 
    nop 
   	nop     
    nop 
    nop 
   	nop     
    nop 
    nop 
    rjmp button_ISR      ; INT2 � ������ ������



INIT:
    ; ������������� ��������� �����
    ldi temp, high(RAMEND)
    out SPH, temp
    ldi temp, low(RAMEND)
    out SPL, temp

    ; ��������� ����� B: LED (PB0) ��� �����, ��������� ��������� � ��� (1)
    ldi temp, 1<<led
    out DDRB, temp
    ldi temp, 0xFF       ; ������������, ��� LED ������ ��� ��������� ���. 1
    out PORTB, temp

    ; ��������� ����� A ��� ������������� ������:
    ; PA0 � SW0, PA1 � SW1, ����� � ����������� ����������� ����������
    clr temp
    out DDRA, temp       ; ���� A ��� ����
    ldi temp, 0b00000011 ; �������� �� PA0 � PA1
    out PORTA, temp

    ; ��������� ����� E ��� INT2 (��������, INT2 �� PE0)
    clr temp
    out DDRE, temp       ; ���� E ��� ����
    ldi temp, (1<<PE0)   ; ��������� ���������� �������� �� PE0
    out PORTE, temp

    ; ���������� ���������� INT2
    ldi temp, (1<<INT2)
    out GICR, temp

    ; ��������� INT2 �� ������������ �� ������� ������ (�������� ��� ��������� �� �����)
    clr temp
    out EMCUCR, temp

    ; ���������� ���������� ����������
    sei

loop:
    rjmp loop

;************************************************************
; ���������� ���������� INT2 (����������� ��� ���� ������)
;************************************************************
button_ISR:
    ; �����������, ����� ������ ������, �� ��������� ������ ����� A.
    ; ���� PA0 = 0 > ������ SW0 ������ (�������� 1 �),
    ; ���� PA0 �� ������, ��������� PA1: ���� PA1 = 0 > SW1 (�������� 2 �).
    in  temp, PINA
    sbic PINA, sw0       ; ���� ��� PA0 ������ (0), ������ SW0 ������
    rjmp check_sw1       ; ���� PA0 �� ������ (��� = 1), ��������� PA1
    rjmp button_sw0      ; ���� PA0 ������, ������� � ��������� SW0

check_sw1:
    sbic PINA, sw1       ; ���� PA1 = 0, ������ SW1 ������
    rjmp wait_release    ; ���� �� ���� ������ �� ������, �������
    rjmp button_sw1      ; ���� PA1 ������, ������� � ��������� SW1

button_sw0:
    ; ��������� ������� ������ SW0: LED ������� �� 1 �
    cbi PORTB, led       ; �������� LED (�������� ��������� � 0)
    rcall delay1
    sbi PORTB, led       ; ��������� LED (������� � ��������� 1)
    rjmp wait_release

button_sw1:
    ; ��������� ������� ������ SW1: LED ������� �� 2 �
    cbi PORTB, led
    rcall delay2
    sbi PORTB, led
    rjmp wait_release

; �������� ���������� ������ (���� PA0 � PA1 �� ������ ����� 1)
wait_release:
wait_rel_loop:
    in temp, PINA
    sbis PINA, sw0       ; ���� PA0 = 1 (��������), ���������� ��������
    rjmp wait_rel_loop   ; ����� � ���
    sbis PINA, sw1
    rjmp wait_rel_loop
    reti                 ; ����� �� ����������

;************************************************************
; ������������ ��������
; delay1 � �������� ~1 �������
; delay2 � �������� ~2 ������� (������� ����� delay1)
;************************************************************
delay1:
    ldi r16, 23
d3:
    ldi r17, 236
d1:
    ldi r18, 245
d2:
    dec r18
    brne d2
    dec r17
    brne d1
    dec r16
    brne d3
    ret

delay2:
    rcall delay1
    rcall delay1
    ret
