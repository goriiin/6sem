.include "m8515def.inc"
.def reg_led = r20
.def temp = r19 
.def store = r21
.def direction = r22

.equ START = 0 ;
.equ STOP = 1 ;
.equ DELAY = 300 ; 300 �� ��������

.org $000
  rjmp INIT

INIT:
  ldi reg_led,0x03   ; ��������� ��������� ���� ���������� ����������� (���� 0 � 1)
  ser temp
  out DDRB,temp      ; ���� B �� �����

  clr temp
  out DDRD,temp      ; ���� D �� ����
  ldi temp,0x03
  out PORTD,temp     ; �������� �� ������

  clr direction      ; ������������� ��������� ����������� �������� (������)

WAITSTART:
  sbic PIND,START    ; ���� ������� ������ START
  rjmp WAITSTART 

LOOP:
  out PORTB,reg_led  ; ����� ��������� �� ���� B

  rcall DELAY_MS     ; �������� 300 ��

  sbic PIND,STOP     ; �������� ������� ������ STOP
  rjmp WAITSTART 

  ; ��������� ����������� ��������
  tst direction      
  breq MOVE_RIGHT   ; ���� direction = 0, ������� ������
  rjmp MOVE_LEFT    ; ���� direction = 1, ������� �����

MOVE_RIGHT:
  lsl reg_led        ; �������� �� 1 ��� �����
  lsl reg_led
  cpi reg_led,0xC0   ; �������� ������� (���� 6 � 7)?
  brne LOOP
  ldi direction,1    ; ������ ����������� �� �����
  rjmp LOOP

MOVE_LEFT:
  lsr reg_led        ; �������� �� 1 ��� ������
  lsr reg_led
  cpi reg_led,0x03   ; �������� ������� (���� 0 � 1)?
  brne LOOP
  clr direction      ; ������ ����������� �� ������
  rjmp LOOP

; ��������� �������� 300 ��
DELAY_MS:
  ldi r16,12
d3: ldi r17,230
d1: ldi r18,245
d2: dec r18
  brne d2
  dec r17
  brne d1
  dec r16
  brne d3
  ret
