.include "m8515def.inc"
.def reg_led = r20      ; ��������� �����������
.def temp = r19
.def direction = r22    ; ����������� ��������: 0-������, 1-�����

.equ START = 0
.equ STOP = 1

.org $000
  rjmp INIT

INIT:
  ldi temp, high(RAMEND)
  out SPH, temp
  ldi temp, low(RAMEND)
  out SPL, temp

  ; ������������� ������
  ldi reg_led, 0b00000011 ; �������� � ���� ���������� ����������� (���� 0 � 1)
  ser temp
  out DDRB, temp          ; ���� B �� �����

  clr temp
  out DDRD, temp          ; ���� D �� ����
  ldi temp, 0x03
  out PORTD, temp         ; �������� �� ������

  ldi direction, 0        ; ��������� ����������� � ������

WAITSTART:
  sbic PIND, START        ; ���� ������� ������ START
  rjmp WAITSTART 

LOOP:
  com reg_led             ; ����������� ��������� ����������� ����� �������!

  out PORTB, reg_led      ; ������� ��������� �� ���� B
  com reg_led             ; ���������� ������� ��� ����������� �������
	
   sbis PIND, STOP         ; �������� ������ STOP
  rjmp WAITSTART 
  rcall DELAY_MS          ; �������� 300 ��

  ; �������� �����������
  tst direction
  breq MOVE_RIGHT
  rjmp MOVE_LEFT

MOVE_RIGHT:
  lsl reg_led             ; �������� ��� ���� �����
 cpi reg_led, 0b11000000 ; ���������, �������� �� ������� (6,7)
  brne LOOP
  ldi direction, 1        ; ������ ����������� �� ����

  sbis PIND, STOP         ; �������� ������ STOP
  rjmp WAITSTART 
  rjmp LOOP

MOVE_LEFT:
  lsr reg_led             ; �������� ��� ���� ������
  cpi reg_led, 0b00000011 ; ���������, �������� �� ������� (0,1)
  brne LOOP
  clr direction           ; ������ ����������� �� �����

   sbis PIND, STOP         ; �������� ������ STOP
  rjmp WAITSTART 
  rjmp LOOP

; ��������� �������� 300 ��
DELAY_MS:
  ldi r16,9
d3: ldi r17,202
d1: ldi r18,221
d2: dec r18
  brne d2
  dec r17
  brne d1
  dec r16
  brne d3
  ret
