   M S B A M L   `   `   `  ����   Wpf3DTube_NonStatic$    Wpf3DTube_NonStatic.MainWindow   \ XPresentationFramework, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35R NWindowsBase, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35W SPresentationCore, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35D 9http://schemas.microsoft.com/winfx/2006/xaml/presentation    R NSystem.Xaml, Version=4.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089:x,http://schemas.microsoft.com/winfx/2006/xaml       � Title$  
MainWindow��$
��1000��$	��525��.��� 
���   � � RowDefinitions �� $��* �.+�'� -    '�Name$ viewport��$��1��$��62,298,-62,-298q� '� Camera =� -   
camera�� � FarPlaneDistance$
 5000Q� � NearPlaneDistance$ 0.125Q� � LookDirection$ -0,-0,-15005� � UpDirection$ 0,1,05� � Position$ 
0,300,40001�	 =� FieldOfView$	 90Q�x�g� -   	model��
 g� Content
 h� -   	group��.����  �� Color$ White��t� $ White�� t� 	Direction$ 0 -1 -15��� -   	spotLight��$ #666666�� �� 	Direction$ 0 -1 -15� �� InnerConeAngle$ 30Q� �� OuterConeAngle$ 60Q� !� Position$ 0 1 301�g� -   FixedShapes��
 h� -   
group2��g� -   MoreFixedShapes��
 h� -	   
group3�� � Triggers :�  :� RoutedEvent$ 
MouseEnter��.[��� .���� -
    g�Name$ myStoryBoardY��