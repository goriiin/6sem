.include "m8515def.inc"
.def reg_led = r20      ; ��������� �����������
.def temp = r19
.def direction = r22    ; ����������� ��������: 0-������, 1-�����

.def temp1 = r16 ;????????? ??????? 1
.def temp2 = r17 ;????????? ??????? 2

.equ START = 0
.equ STOP = 1

.macro reload_timer ;??????
	 ldi temp1,high(48239) ; ??? ???????? ?
	 ldi temp2,low(48239) ; ??????
	 out TCNT1H,temp1 ; ????????? ????????
	 out TCNT1L,temp2 ; ??? ?????? ???????
.endmacro

.org $000
  	rjmp INIT

.org $0006
	 rjmp T1_OVF ;

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

	ldi temp1,(1<<TOIE1) ;?????????? ??????????
	out TIMSK,temp1 ; ??????? ?? ????????????

	ser temp
	out PIND, temp

	reload_timer
	ldi temp1, ((1<<CS11)|(1<<CS10))  ; ����� �������� 64
	out TCCR1B, temp1
	
	sei
	

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
	clt

WAIT_END:
	brtc WAIT_END

	ret

; ���������� ���������� ������������ �������
T1_OVF:

	reload_timer
	ldi temp1, ((1<<CS11)|(1<<CS10))  ; ����� �������� 64
	out TCCR1B, temp1

	set

	reti
